import requests
import os
import json
import numpy as np
import pandas as pd
from .config import (CACHE_DIR, TARGET_TISSUES, GTEX_UTERUS_TISSUE_ID,
                     ENABLE_CELLXGENE, CELLXGENE_DATA_PATH, CELLXGENE_REPRODUCTIVE_TISSUES)


class ExternalDataManager:
    """
    V2: Handles GTEx and CellxGene external data with local caching.

    GTEx off-target metric: max log2((TPM_off_target + 1) / (TPM_uterus + 1))
    across non-reproductive tissues. Uterus is used as reference so the metric is
    comparable in scale to EuE-normalized logFC values from Tan et al. 2022.
    Note: GTEx 'Uterus' is bulk mixed tissue (myometrium + endometrium) and is not
    equivalent to EuE (scRNA-seq eutopic endometrial cells). Values are directionally
    valid but not directly interchangeable with Supp Table 5 logFC.

    CellxGene off-target metric: max % cells expressing the gene across non-reproductive
    healthy tissues. Measures breadth of off-target expression. Loaded from a
    user-exported CSV (cellxgene_data.csv).
    """

    GTEX_API_BASE = "https://gtexportal.org/api/v2"

    def __init__(self):
        os.makedirs(CACHE_DIR, exist_ok=True)
        self.gtex_cache_path  = os.path.join(CACHE_DIR, "gtex_cache_v2.json")
        self.gtex_id_map_path = os.path.join(CACHE_DIR, "gtex_id_map.json")
        self.gtex_cache  = self._load_cache(self.gtex_cache_path)
        self.gtex_id_map = self._load_cache(self.gtex_id_map_path)

        # CellxGene: load from user-exported CSV
        self.cellxgene_df = None
        if ENABLE_CELLXGENE:
            if os.path.exists(CELLXGENE_DATA_PATH):
                self.cellxgene_df = pd.read_csv(CELLXGENE_DATA_PATH)
                print(f"  CellxGene data loaded: {len(self.cellxgene_df)} rows")
            else:
                print(f"  Warning: ENABLE_CELLXGENE=True but '{CELLXGENE_DATA_PATH}' not found. Skipping.")

    # ------------------------------------------------------------------ #
    #  Cache helpers
    # ------------------------------------------------------------------ #
    def _load_cache(self, path):
        if os.path.exists(path):
            with open(path, 'r') as f:
                return json.load(f)
        return {}

    def _save_cache(self, cache, path):
        with open(path, 'w') as f:
            json.dump(cache, f)

    # ------------------------------------------------------------------ #
    #  GTEx
    # ------------------------------------------------------------------ #
    def get_gencode_id(self, symbol):
        """Maps gene symbol to GTEx Gencode ID."""
        if symbol in self.gtex_id_map:
            return self.gtex_id_map[symbol]
        try:
            r = requests.get(
                f"{self.GTEX_API_BASE}/reference/gene",
                params={"geneId": symbol, "format": "json"}, timeout=10)
            if r.status_code == 200:
                data = r.json().get('data', [])
                if data:
                    gid = data[0].get('gencodeId')
                    self.gtex_id_map[symbol] = gid
                    self._save_cache(self.gtex_id_map, self.gtex_id_map_path)
                    return gid
        except Exception as e:
            print(f"    ! GTEx ID lookup failed for {symbol}: {e}")
        return None

    def _fetch_gtex_tissue_tpm(self, gene_symbol):
        """Returns {tissue_id: median_tpm} for all GTEx tissues. Cached."""
        cache_key = f"v2_{gene_symbol}"
        if cache_key in self.gtex_cache:
            return self.gtex_cache[cache_key]

        gencode_id = self.get_gencode_id(gene_symbol)
        if not gencode_id:
            return {}
        try:
            r = requests.get(
                f"{self.GTEX_API_BASE}/expression/medianGeneExpression",
                params={"datasetId": "gtex_v8", "gencodeId": gencode_id, "format": "json"},
                timeout=15)
            if r.status_code == 200:
                tissue_tpm = {
                    e.get('tissueSiteDetailId', ''): e.get('median', 0.0)
                    for e in r.json().get('data', [])
                }
                self.gtex_cache[cache_key] = tissue_tpm
                self._save_cache(self.gtex_cache, self.gtex_cache_path)
                return tissue_tpm
        except Exception as e:
            print(f"    ! GTEx expression fetch failed for {gene_symbol}: {e}")
        return {}

    def get_off_target_burden(self, gene_symbol):
        """
        V2 GTEx metric: max log2((TPM_off_target + 1) / (TPM_uterus + 1))
        across non-reproductive tissues. Clipped at 0 — a negative ratio means
        the gene is less expressed there than in uterus, which is fine.
        Returns np.nan if no GTEx data is available.
        """
        tissue_tpm = self._fetch_gtex_tissue_tpm(gene_symbol)
        if not tissue_tpm:
            return np.nan

        # Get uterus reference TPM
        uterus_tpm = 0.0
        for tissue_id, tpm in tissue_tpm.items():
            if GTEX_UTERUS_TISSUE_ID.lower() in tissue_id.lower():
                uterus_tpm = tpm
                break

        ratios = []
        for tissue_id, tpm in tissue_tpm.items():
            is_reproductive = any(
                tt.lower().replace(' ', '_') in tissue_id.lower()
                for tt in TARGET_TISSUES
            )
            if not is_reproductive:
                ratio = np.log2((tpm + 1) / (uterus_tpm + 1))
                ratios.append(max(0.0, ratio))  # Clip at 0: negative = less than uterus = fine

        return max(ratios) if ratios else np.nan

    # ------------------------------------------------------------------ #
    #  CellxGene
    # ------------------------------------------------------------------ #
    def get_cellxgene_burden(self, gene_symbol):
        """
        Returns max % of cells expressing the gene across non-reproductive
        healthy tissues. Reads from user-exported cellxgene_data.csv.
        Expected columns: Gene, Tissue, Percent_Cells
        """
        if self.cellxgene_df is None:
            return np.nan

        gene_rows = self.cellxgene_df[
            self.cellxgene_df['Gene'].str.upper() == gene_symbol.upper()
        ]
        if gene_rows.empty:
            return np.nan

        non_repro = gene_rows[
            ~gene_rows['Tissue'].str.lower().isin(CELLXGENE_REPRODUCTIVE_TISSUES)
        ]
        return float(non_repro['Percent_Cells'].max()) if not non_repro.empty else 0.0

    # ------------------------------------------------------------------ #
    #  Batch
    # ------------------------------------------------------------------ #
    def fetch_batch_specificity(self, gene_list):
        """Returns {gene: {'gtex_burden': float, 'cellxgene_burden': float}}"""
        results = {}
        total = len(gene_list)
        for i, gene in enumerate(gene_list, 1):
            print(f"  [{i}/{total}] {gene}")
            results[gene] = {
                'gtex_burden':      self.get_off_target_burden(gene),
                'cellxgene_burden': self.get_cellxgene_burden(gene)
            }
        return results
