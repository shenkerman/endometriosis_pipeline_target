import pandas as pd
import numpy as np
from .config import (CELL_TYPES_INFO, ENABLE_EXTERNAL_SPECIFICITY, TOP_N_FOR_EXTERNAL,
                     CPM_PERCENTILE_THRESHOLD, EUE_NEAR_ZERO_THRESHOLD, ENABLE_CELLXGENE)
from .external_data import ExternalDataManager


class PipelineDataProcessor:
    def __init__(self, file_path):
        self.file_path = file_path
        self.calculated_weights = {}
        self.external_manager = ExternalDataManager() if ENABLE_EXTERNAL_SPECIFICITY else None

    # ------------------------------------------------------------------ #
    #  Data Loading
    # ------------------------------------------------------------------ #
    def load_raw_data(self):
        """Loads the Excel file and parses it into a single combined DataFrame."""
        df_raw = pd.read_excel(self.file_path, header=None)
        header_row = 5
        data_start = 6
        all_dfs = []

        for start_col, name in CELL_TYPES_INFO:
            end_col = start_col + 9
            cell_df = df_raw.iloc[data_start:, start_col:end_col + 1].copy()
            cols = df_raw.iloc[header_row, start_col:end_col + 1].tolist()
            cell_df.columns = cols
            cell_df = cell_df[cell_df['Gene'].notna()]
            cell_df['Cell Type'] = name
            all_dfs.append(cell_df)

        return pd.concat(all_dfs, ignore_index=True)

    # ------------------------------------------------------------------ #
    #  Filters — F1 and F2 only (F4 removed in V2)
    # ------------------------------------------------------------------ #
    @staticmethod
    def apply_filters(df):
        """
        Applies F1 (logCPM) and F2 (FDR) filters.

        F4 'Inconsistent direction' has been removed in V2.
        A gene downregulated in EuE but upregulated in lesions is a desirable
        biomarker profile, not a reason to drop the gene.
        """
        df = df.copy()
        df['Status'] = 'PASS'
        df['Dropped Reason'] = ''

        def filter_row(row):
            reasons = []
            if pd.isna(row['logCPM']) or row['logCPM'] < 1:
                reasons.append('F1 (logCPM < 1)')
            if pd.isna(row['FDR']) or row['FDR'] > 0.01:
                reasons.append('F2 (FDR > 0.01)')
            return ('DROPPED', '; '.join(reasons)) if reasons else ('PASS', '')

        results = df.apply(filter_row, axis=1)
        df['Status']         = [r[0] for r in results]
        df['Dropped Reason'] = [r[1] for r in results]
        return df

    # ------------------------------------------------------------------ #
    #  Entropy Weighting
    # ------------------------------------------------------------------ #
    @staticmethod
    def _calculate_entropy_weights(matrix):
        """Calculates weights using the Entropy Weighting Method."""
        n = matrix.shape[0]
        # Entropy is undefined for a single row (log(1)=0); fall back to equal weights
        if n <= 1:
            return np.ones(matrix.shape[1]) / matrix.shape[1]
        matrix = np.abs(matrix) + 1e-9
        p = matrix / matrix.sum(axis=0)
        k = 1.0 / np.log(n)
        entropy = -k * np.nan_to_num(p * np.log(p)).sum(axis=0)
        d = 1.0 - entropy
        return np.ones(matrix.shape[1]) / matrix.shape[1] if d.sum() == 0 else d / d.sum()

    # ------------------------------------------------------------------ #
    #  TOPSIS
    # ------------------------------------------------------------------ #
    def _run_topsis(self, df, criteria_cols, weights, minimize_cols=None):
        """
        Runs TOPSIS on the given criteria columns.
        minimize_cols: list of column names where a lower value is better
                       (e.g. off-target burden). All other columns are maximized.
        """
        if minimize_cols is None:
            minimize_cols = []

        matrix = df[criteria_cols].apply(pd.to_numeric, errors='coerce').fillna(0).values
        norms = np.sqrt((matrix ** 2).sum(axis=0)) + 1e-9
        norm_matrix = matrix / norms
        weighted = norm_matrix * weights

        v_ideal = np.zeros(len(criteria_cols))
        v_anti  = np.zeros(len(criteria_cols))
        for i, col in enumerate(criteria_cols):
            if col in minimize_cols:
                v_ideal[i] = weighted[:, i].min()   # Best = smallest
                v_anti[i]  = weighted[:, i].max()   # Worst = largest
            else:
                v_ideal[i] = weighted[:, i].max()   # Best = largest
                v_anti[i]  = weighted[:, i].min()   # Worst = smallest

        s_ideal = np.sqrt(((weighted - v_ideal) ** 2).sum(axis=1))
        s_anti  = np.sqrt(((weighted - v_anti)  ** 2).sum(axis=1))
        return s_anti / (s_ideal + s_anti + 1e-9)

    # ------------------------------------------------------------------ #
    #  Derived Columns
    # ------------------------------------------------------------------ #
    @staticmethod
    def _add_derived_columns(df):
        """
        Adds:
        - CPM: back-transformed from logCPM (2^logCPM), rounded to 2 decimal places
        - CPM_percentile: rank percentile within the passed gene set (0-1)
        - EuE_near_zero: True if |EuE logFC| < EUE_NEAR_ZERO_THRESHOLD (flag for manual review)
        """
        df = df.copy()

        # CPM back-transform
        df['CPM'] = (2 ** pd.to_numeric(df['logCPM'], errors='coerce')).round(2)

        # CPM percentile within the current result set
        df['CPM_percentile'] = pd.to_numeric(df['CPM'], errors='coerce').rank(pct=True).round(3)

        # EuE near-zero flag
        eue_col = 'logFC.eueVSctrl'
        if eue_col in df.columns:
            df['EuE_near_zero'] = pd.to_numeric(df[eue_col], errors='coerce').abs() < EUE_NEAR_ZERO_THRESHOLD
        else:
            df['EuE_near_zero'] = np.nan

        return df

    # ------------------------------------------------------------------ #
    #  Ranking
    # ------------------------------------------------------------------ #
    def calculate_ranking(self, df):
        """
        V2 two-stage ranking:

        Stage 1 (Local): TOPSIS on lesion-specificity distances:
            specificity_ecp = |EcP - EuE|
            specificity_eco = |EcO - EuE|
        EcPA is excluded. Direction (up/down) does not affect ranking —
        both are valid biomarker profiles.

        Stage 2 (External): Refines top TOP_N_FOR_EXTERNAL genes with:
            gtex_burden      = max log2((TPM_off_target+1)/(TPM_uterus+1))
            cellxgene_burden = max % cells expressing across non-reproductive tissues
        Both are minimize criteria. No penalty for genes outside top N.
        """
        df = df.copy()
        passed  = df[df['Status'] == 'PASS'].copy()
        dropped = df[df['Status'] == 'DROPPED'].copy()

        if passed.empty:
            self.calculated_weights = {}
            return df

        # -- Compute specificity distances (V2 core change) --
        eue = pd.to_numeric(passed.get('logFC.eueVSctrl', pd.Series(np.nan, index=passed.index)), errors='coerce')
        ecp = pd.to_numeric(passed.get('logFC.ecpVSctrl', pd.Series(np.nan, index=passed.index)), errors='coerce')
        eco = pd.to_numeric(passed.get('logFC.ecoVSctrl', pd.Series(np.nan, index=passed.index)), errors='coerce')

        passed['specificity_ecp'] = (ecp - eue).abs()
        passed['specificity_eco'] = (eco - eue).abs()
        # EcPA intentionally excluded from all criteria

        # -- Stage 1: Local ranking --
        print("Stage 1: Ranking on lesion-specificity distances (|EcP-EuE|, |EcO-EuE|)...")
        local_criteria = ['specificity_ecp', 'specificity_eco']
        local_matrix   = passed[local_criteria].apply(pd.to_numeric).fillna(0).values
        local_weights  = self._calculate_entropy_weights(local_matrix)
        passed['Local_Score'] = self._run_topsis(passed, local_criteria, local_weights)
        self.calculated_weights = {'Stage 1 (Local)': dict(zip(local_criteria, local_weights))}

        # -- Stage 2: External validation --
        if ENABLE_EXTERNAL_SPECIFICITY:
            print(f"Stage 2: External validation for top {TOP_N_FOR_EXTERNAL} genes...")
            top_genes    = passed.sort_values('Local_Score', ascending=False)['Gene'].unique()[:TOP_N_FOR_EXTERNAL]
            external_raw = self.external_manager.fetch_batch_specificity(top_genes)

            passed['gtex_burden']      = passed['Gene'].map({g: v['gtex_burden']      for g, v in external_raw.items()})
            passed['cellxgene_burden'] = passed['Gene'].map({g: v['cellxgene_burden'] for g, v in external_raw.items()})

            refined = passed[passed['gtex_burden'].notna()].copy()
            other   = passed[passed['gtex_burden'].isna()].copy()

            if not refined.empty:
                off_target_cols = ['gtex_burden']
                if ENABLE_CELLXGENE and refined['cellxgene_burden'].notna().any():
                    off_target_cols.append('cellxgene_burden')

                full_criteria   = local_criteria + off_target_cols
                full_matrix     = refined[full_criteria].apply(pd.to_numeric).fillna(0).values
                refined_weights = self._calculate_entropy_weights(full_matrix)
                self.calculated_weights['Stage 2 (Refined)'] = dict(zip(full_criteria, refined_weights))

                refined['Score'] = self._run_topsis(
                    refined, full_criteria, refined_weights,
                    minimize_cols=off_target_cols
                )

                # off_target_agree: True if GTEx and CellxGene agree on risk level
                if 'cellxgene_burden' in off_target_cols:
                    gtex_med = refined['gtex_burden'].median()
                    cxg_med  = refined['cellxgene_burden'].median()
                    refined['off_target_agree'] = (
                        ((refined['gtex_burden'] < gtex_med) & (refined['cellxgene_burden'] < cxg_med)) |
                        ((refined['gtex_burden'] >= gtex_med) & (refined['cellxgene_burden'] >= cxg_med))
                    )

                # V2: no penalty for genes outside top N — they keep their local score
                other['Score'] = other['Local_Score']

                passed = pd.concat([refined, other])
            else:
                passed['Score'] = passed['Local_Score']
        else:
            passed['Score'] = passed['Local_Score']

        # -- Add derived columns (CPM, percentile, EuE flag) --
        passed = self._add_derived_columns(passed)

        # -- Backward compat alias for visualizer --
        if 'gtex_burden' in passed.columns:
            passed['off_target_burden'] = passed['gtex_burden']

        # -- Final sort: top 15% CPM first, then by Score within each group --
        high_cpm    = passed['CPM_percentile'] >= CPM_PERCENTILE_THRESHOLD
        passed_high = passed[high_cpm].sort_values('Score', ascending=False)
        passed_low  = passed[~high_cpm].sort_values('Score', ascending=False)
        passed = pd.concat([passed_high, passed_low])
        passed['Rank'] = range(1, len(passed) + 1)

        dropped['Rank']        = np.nan
        dropped['Score']       = np.nan
        dropped['Local_Score'] = np.nan

        return pd.concat([passed, dropped], ignore_index=True)
