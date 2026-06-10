import requests
import os
import json
import pickle
import numpy as np
import pandas as pd
from .config import (CACHE_DIR, TARGET_TISSUES, GTEX_UTERUS_TISSUE_ID,
                     ENABLE_CELLXGENE, CELLXGENE_REPRODUCTIVE_TISSUES)


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
    healthy tissues (disease == 'normal'). Fetched via cellxgene-census Python API.
    Results are cached locally in .api_cache/cellxgene_cache.pkl.
    """

    GTEX_API_BASE = "https://gtexportal.org/api/v2"

    def __init__(self):
        os.makedirs(CACHE_DIR, exist_ok=True)
        self.gtex_cache_path       = os.path.join(CACHE_DIR, "gtex_cache_v2.json")
        self.gtex_id_map_path      = os.path.join(CACHE_DIR, "gtex_id_map.json")
        self.cellxgene_cache_path  = os.path.join(CACHE_DIR, "cellxgene_cache.pkl")

        self.gtex_cache       = self._load_json(self.gtex_cache_path)
        self.gtex_id_map      = self._load_json(self.gtex_id_map_path)
        self.cellxgene_cache  = self._load_pickle(self.cellxgene_cache_path)

        if ENABLE_CELLXGENE:
            try:
                import cellxgene_census  # noqa: F401
            except ImportError:
                print("  Warning: cellxgene-census not installed. "
                      "Run: pip install cellxgene-census")
                print("  CellxGene off-target will be skipped.")

    # ------------------------------------------------------------------ #
    #  Cache helpers
    # ------------------------------------------------------------------ #
    def _load_json(self, path):
        if os.path.exists(path):
            with open(path, 'r') as f:
                return json.load(f)
        return {}

    def _save_json(self, data, path):
        with open(path, 'w') as f:
            json.dump(data, f)

    def _load_pickle(self, path):
        if os.path.exists(path):
            with open(path, 'rb') as f:
                return pickle.load(f)
        return {}

    def _save_pickle(self, data, path):
        with open(path, 'wb') as f:
            pickle.dump(data, f)

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
                    self._save_json(self.gtex_id_map, self.gtex_id_map_path)
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
                self._save_json(self.gtex_cache, self.gtex_cache_path)
                return tissue_tpm
        except Exception as e:
            print(f"    ! GTEx expression fetch failed for {gene_symbol}: {e}")
        return {}

    def get_off_target_burden(self, gene_symbol):
        """
        V2 GTEx metric: max log2((TPM_off_target + 1) / (TPM_uterus + 1))
        across non-reproductive tissues. Clipped at 0 - a negative ratio means
        the gene is less expressed there than in uterus, which is fine.
        Returns np.nan if no GTEx data is available.
        """
        tissue_tpm = self._fetch_gtex_tissue_tpm(gene_symbol)
        if not tissue_tpm:
            return np.nan

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
                ratios.append(max(0.0, ratio))

        return max(ratios) if ratios else np.nan

    # ------------------------------------------------------------------ #
    #  CellxGene Census API
    # ------------------------------------------------------------------ #
    def _fetch_cellxgene_batch_api(self, gene_list):
        """
        Query the CellxGene Census for % cells expressing each gene
        across healthy (disease == 'normal') non-reproductive tissues.
        Results stored in self.cellxgene_cache and saved to disk.
        First call per gene set may take several minutes; subsequent runs use cache.
        """
        try:
            import cellxgene_census
            import scipy.sparse as sp
        except ImportError:
            print("  ! cellxgene-census not installed - skipping CellxGene burden.")
            return

        to_fetch = [g for g in gene_list if g not in self.cellxgene_cache]
        if not to_fetch:
            return

        print(f"  Fetching CellxGene Census data for {len(to_fetch)} genes "
              f"(this may take a few minutes the first time)...")

        try:
            var_filter = " or ".join([f"feature_name == '{g}'" for g in to_fetch])

            with cellxgene_census.open_soma() as census:
                adata = cellxgene_census.get_anndata(
                    census,
                    organism="Homo sapiens",
                    obs_value_filter="disease == 'normal' and is_primary_data == True",
                    var_value_filter=var_filter,
                    obs_column_names=["tissue_general"],
                )

            if adata.n_vars == 0 or adata.n_obs == 0:
                print("  ! CellxGene: no data returned for query.")
                return

            tissues = adata.obs['tissue_general'].values
            unique_tissues = np.unique(tissues)
            X = adata.X

            for i, gene in enumerate(adata.var['feature_name'].values):
                if sp.issparse(X):
                    col = np.asarray(X[:, i].todense()).flatten()
                else:
                    col = np.asarray(X[:, i]).flatten()

                expressing = col > 0

                non_repro_pcts = []
                for t in unique_tissues:
                    if any(r in t.lower() for r in CELLXGENE_REPRODUCTIVE_TISSUES):
                        continue
                    mask = tissues == t
                    pct = float(expressing[mask].mean()) if mask.any() else 0.0
                    non_repro_pcts.append(pct)

                self.cellxgene_cache[gene] = (
                    float(max(non_repro_pcts)) if non_repro_pcts else np.nan
                )

            self._save_pickle(self.cellxgene_cache, self.cellxgene_cache_path)
            print(f"  CellxGene data fetched and cached for {len(to_fetch)} genes.")

        except Exception as e:
            print(f"  ! CellxGene Census fetch error: {e}")

    def get_cellxgene_burden(self, gene_symbol):
        """Returns cached % expressing value for a single gene."""
        if not ENABLE_CELLXGENE:
            return np.nan
        return self.cellxgene_cache.get(gene_symbol, np.nan)

    # ------------------------------------------------------------------ #
    #  Batch entry point (called by processor)
    # ------------------------------------------------------------------ #
    def fetch_batch_specificity(self, gene_list):
        """
        Returns {gene: {'gtex_burden': float, 'cellxgene_burden': float}}
        GTEx: fetched gene-by-gene with JSON cache.
        CellxGene: fetched in one batch via Census API with pickle cache.
        """
        # Pre-fetch CellxGene for all genes at once (single Census connection)
        if ENABLE_CELLXGENE:
            self._fetch_cellxgene_batch_api(gene_list)

        results = {}
        total = len(gene_list)
        for i, gene in enumerate(gene_list, 1):
            print(f"  [{i}/{total}] GTEx: {gene}")
            results[gene] = {
                'gtex_burden':      self.get_off_target_burden(gene),
                'cellxgene_burden': self.get_cellxgene_burden(gene),
            }
        return results
