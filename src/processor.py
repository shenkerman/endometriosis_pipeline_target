import pandas as pd
import numpy as np
from .config import (CELL_TYPES_INFO, ENABLE_EXTERNAL_SPECIFICITY, TOP_N_FOR_EXTERNAL,
                     CPM_PERCENTILE_THRESHOLD, ENABLE_CELLXGENE)
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
        - CPM_percentile: rank percentile within each cell type (not cross-cell-type)
        """
        df = df.copy()

        # CPM back-transform
        df['CPM'] = (2 ** pd.to_numeric(df['logCPM'], errors='coerce')).round(2)

        # CPM percentile within each cell type (not cross-cell-type)
        # Different cell types have very different baseline CPMs (e.g. cDC2 ~149 vs Treg ~367),
        # so cross-cell-type ranking would be misleading as a tiebreaker.
        if 'Cell Type' in df.columns:
            df['CPM_percentile'] = (
                df.groupby('Cell Type')['CPM']
                .transform(lambda x: pd.to_numeric(x, errors='coerce').rank(pct=True))
                .round(3)
            )
        else:
            df['CPM_percentile'] = pd.to_numeric(df['CPM'], errors='coerce').rank(pct=True).round(3)

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

        Stage 2 (External, info-only): Fetches off-target data for top TOP_N_FOR_EXTERNAL genes:
            gtex_burden_vs_uterus = max log2((TPM_off_target+1)/(TPM_uterus+1))
            cellxgene_burden      = max % cells expressing across non-reproductive tissues
        These are stored as output columns for manual review but do NOT affect Score.
        Final Score = Local_Score for all genes.
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

        # -- Stage 2: External data — informational only, does NOT affect ranking --
        # GTEx and CellxGene burden columns are fetched for the top N genes and stored
        # as output columns for manual review. Score remains Local_Score for all genes.
        # Rationale: off-target data is on a different scale from lesion logFC and
        # introduces cross-source normalization assumptions that should be a human
        # decision, not an automated TOPSIS weight.
        if ENABLE_EXTERNAL_SPECIFICITY:
            print(f"Stage 2: Fetching external data for top {TOP_N_FOR_EXTERNAL} genes "
                  f"(info only — does not change ranking)...")
            top_genes    = passed.sort_values('Local_Score', ascending=False)['Gene'].unique()[:TOP_N_FOR_EXTERNAL]
            external_raw = self.external_manager.fetch_batch_specificity(top_genes)

            passed['gtex_burden_vs_uterus'] = passed['Gene'].map(
                {g: v['gtex_burden_vs_uterus'] for g, v in external_raw.items()})
            passed['cellxgene_burden'] = passed['Gene'].map(
                {g: v['cellxgene_burden'] for g, v in external_raw.items()})

            # off_target_agree: informational flag — True if GTEx and CellxGene agree on risk
            has_gtex = passed['gtex_burden_vs_uterus'].notna().any()
            has_cxg  = passed['cellxgene_burden'].notna().any()
            if has_gtex and has_cxg:
                gtex_med = passed['gtex_burden_vs_uterus'].median()
                cxg_med  = passed['cellxgene_burden'].median()
                passed['off_target_agree'] = (
                    ((passed['gtex_burden_vs_uterus'] < gtex_med) & (passed['cellxgene_burden'] < cxg_med)) |
                    ((passed['gtex_burden_vs_uterus'] >= gtex_med) & (passed['cellxgene_burden'] >= cxg_med))
                )

        # Ranking is based solely on Stage 1 lesion-specificity score
        passed['Score'] = passed['Local_Score']

        # -- Add derived columns (CPM, percentile, EuE flag) --
        passed = self._add_derived_columns(passed)

        # -- Backward compat alias for visualizer --
        if 'gtex_burden_vs_uterus' in passed.columns:
            passed['off_target_burden'] = passed['gtex_burden_vs_uterus']

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
