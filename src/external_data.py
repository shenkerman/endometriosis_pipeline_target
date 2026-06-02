import requests
import os
import json
import pandas as pd
from .config import CACHE_DIR, TARGET_TISSUES

class ExternalDataManager:
    """Handles interaction with external bioinformatics APIs with local caching."""
    
    GTEx_API_BASE = "https://gtexportal.org/api/v2"

    def __init__(self):
        os.makedirs(CACHE_DIR, exist_ok=True)
        self.gtex_cache_path = os.path.join(CACHE_DIR, "gtex_cache.json")
        self.gtex_id_map_path = os.path.join(CACHE_DIR, "gtex_id_map.json")
        self.gtex_cache = self._load_cache(self.gtex_cache_path)
        self.gtex_id_map = self._load_cache(self.gtex_id_map_path)

    def _load_cache(self, path):
        if os.path.exists(path):
            with open(path, 'r') as f:
                return json.load(f)
        return {}

    def _save_cache(self, cache, path):
        with open(path, 'w') as f:
            json.dump(cache, f)

    def get_gencode_id(self, symbol):
        """Maps Gene Symbol to GTEx-compatible Gencode ID (Ensembl ID)."""
        if symbol in self.gtex_id_map:
            return self.gtex_id_map[symbol]
        
        print(f"    > Searching GTEx for Gencode ID: {symbol}...")
        url = f"{self.GTEx_API_BASE}/reference/gene"
        params = {"geneId": symbol, "format": "json"}
        
        try:
            response = requests.get(url, params=params, timeout=10)
            if response.status_code == 200:
                data = response.json().get('data', [])
                if data:
                    gencode_id = data[0].get('gencodeId')
                    self.gtex_id_map[symbol] = gencode_id
                    self._save_cache(self.gtex_id_map, self.gtex_id_map_path)
                    return gencode_id
        except Exception as e:
            print(f"      ! Error mapping symbol {symbol}: {e}")
        return None

    def get_off_target_burden(self, gene_symbol):
        """
        Fetches median expression from GTEx and calculates total burden 
        outside of target reproductive tissues.
        """
        if gene_symbol in self.gtex_cache:
            return self.gtex_cache[gene_symbol]

        gencode_id = self.get_gencode_id(gene_symbol)
        if not gencode_id:
            return 0.0

        print(f"  > Fetching GTEx expression data for {gene_symbol} ({gencode_id})...")
        url = f"{self.GTEx_API_BASE}/expression/medianGeneExpression"
        params = {
            "datasetId": "gtex_v8",
            "gencodeId": gencode_id,
            "format": "json"
        }
        
        try:
            response = requests.get(url, params=params, timeout=15)
            if response.status_code == 200:
                data = response.json().get('data', [])
                if not data:
                    return 0.0
                
                burden = 0.0
                for entry in data:
                    tissue = entry.get('tissueSiteDetailId', '')
                    # Target tissues are Ovary, Uterus, Vagina, Fallopian Tube, Cervix
                    is_target = any(tt.lower().replace(' ', '_') in tissue.lower() for tt in TARGET_TISSUES)
                    if not is_target:
                        burden += entry.get('median', 0.0)
                
                self.gtex_cache[gene_symbol] = burden
                self._save_cache(self.gtex_cache, self.gtex_cache_path)
                return burden
        except Exception as e:
            print(f"    ! Error fetching GTEx for {gene_symbol}: {e}")
            
        return 0.0

    def fetch_batch_specificity(self, gene_list):
        """Enriches a gene list with external specificity data."""
        results = {}
        for gene in gene_list:
            results[gene] = self.get_off_target_burden(gene)
        return results
