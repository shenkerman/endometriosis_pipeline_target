import pandas as pd
import numpy as np
from .config import CELL_TYPES_INFO, RANKING_WEIGHTS, ENABLE_EXTERNAL_SPECIFICITY, TOP_N_FOR_EXTERNAL
from .external_data import ExternalDataManager

class PipelineDataProcessor:
    def __init__(self, file_path):
        self.file_path = file_path
        self.calculated_weights = {}
        self.external_manager = ExternalDataManager() if ENABLE_EXTERNAL_SPECIFICITY else None

    def load_raw_data(self):
        """Loads the Excel file and parses it into a single combined DataFrame."""
        df_raw = pd.read_excel(self.file_path, header=None)
        header_row = 5
        data_start_row = 6
        
        all_dfs = []
        for start_col, name in CELL_TYPES_INFO:
            end_col = start_col + 9
            cell_df = df_raw.iloc[data_start_row:, start_col:end_col+1].copy()
            
            # Set columns from row 5
            cols = df_raw.iloc[header_row, start_col:end_col+1].tolist()
            cell_df.columns = cols
            
            # Clean up and annotate
            cell_df = cell_df[cell_df['Gene'].notna()]
            cell_df['Cell Type'] = name
            all_dfs.append(cell_df)
            
        return pd.concat(all_dfs, ignore_index=True)

    @staticmethod
    def apply_filters(df):
        """Applies filters F1, F2, F4 and sets the Status."""
        df = df.copy()
        df['Status'] = 'PASS'
        df['Dropped Reason'] = ''
        
        logfc_cols = [c for c in df.columns if isinstance(c, str) and 'logFC' in c]
        
        def filter_row(row):
            reasons = []
            if pd.isna(row['logCPM']) or row['logCPM'] < 1:
                reasons.append('F1 (logCPM < 1)')
            if pd.isna(row['FDR']) or row['FDR'] > 0.01:
                reasons.append('F2 (FDR > 0.01)')
            
            vals = [row[c] for c in logfc_cols if pd.notna(row[c])]
            if not vals:
                reasons.append('F4 (No logFC data)')
            elif not (all(v > 0 for v in vals) or all(v < 0 for v in vals)):
                reasons.append('F4 (Inconsistent direction)')
            
            if reasons:
                return 'DROPPED', '; '.join(reasons)
            return 'PASS', ''

        status_reasons = df.apply(filter_row, axis=1)
        df['Status'] = [x[0] for x in status_reasons]
        df['Dropped Reason'] = [x[1] for x in status_reasons]
        return df

    @staticmethod
    def _calculate_entropy_weights(matrix):
        """Calculates weights using the Entropy Weighting Method."""
        # Ensure positive values and avoid log(0)
        matrix = np.abs(matrix) + 1e-9
        column_sums = matrix.sum(axis=0)
        p_matrix = matrix / column_sums
        
        n = matrix.shape[0]
        k = 1.0 / np.log(n)
        
        p_ln_p = p_matrix * np.log(p_matrix)
        p_ln_p = np.nan_to_num(p_ln_p)
        entropy = -k * p_ln_p.sum(axis=0)
        
        diversification = 1.0 - entropy
        # In case diversification is 0 for a column (all values same)
        if diversification.sum() == 0:
            return np.ones(matrix.shape[1]) / matrix.shape[1]
            
        weights = diversification / diversification.sum()
        return weights

    def _run_topsis(self, df, criteria_cols, weights, specificity_cols=[]):
        """Internal helper to run TOPSIS on a specific matrix."""
        # 1. Normalize
        matrix = df[criteria_cols].apply(pd.to_numeric, errors='coerce').fillna(0).abs().values
        sums_of_squares = (matrix**2).sum(axis=0)
        norm_matrix = matrix / np.sqrt(sums_of_squares + 1e-9)
        
        # 2. Weighted Normalized Matrix
        weighted_matrix = norm_matrix * weights
        
        # 3. Ideal and Anti-Ideal
        v_ideal = []
        v_anti_ideal = []
        for i, col in enumerate(criteria_cols):
            if col in ['logFC.eueVSctrl', 'off_target_burden'] or col in specificity_cols:
                # Minimize (Ideal is 0)
                v_ideal.append(weighted_matrix[:, i].min())
                v_anti_ideal.append(weighted_matrix[:, i].max())
            else:
                # Maximize (Ideal is max logFC)
                v_ideal.append(weighted_matrix[:, i].max())
                v_anti_ideal.append(weighted_matrix[:, i].min())
        
        v_ideal = np.array(v_ideal)
        v_anti_ideal = np.array(v_anti_ideal)
        
        # 4. Geometric Distance
        s_ideal = np.sqrt(((weighted_matrix - v_ideal)**2).sum(axis=1))
        s_anti_ideal = np.sqrt(((weighted_matrix - v_anti_ideal)**2).sum(axis=1))
        
        # 5. Score
        return s_anti_ideal / (s_ideal + s_anti_ideal)

    def calculate_ranking(self, df):
        """
        Calculates ranking in two stages:
        Stage 1: Identify interesting genes based on local dataset.
        Stage 2: Refine ranking for top candidates using GTEx global specificity.
        """
        df = df.copy()
        passed = df[df['Status'] == 'PASS'].copy()
        dropped = df[df['Status'] == 'DROPPED'].copy()
        
        if passed.empty:
            self.calculated_weights = {}
            return df

        # --- Stage 1: Local Dataset Ranking ---
        local_criteria = ['logFC.eueVSctrl', 'logFC.ecpVSctrl', 'logFC.ecoVSctrl', 'logFC.ecpaVSctrl']
        local_criteria = [c for c in local_criteria if c in passed.columns]
        
        print("Stage 1: Ranking based on local study data...")
        local_matrix = passed[local_criteria].apply(pd.to_numeric).fillna(0).abs().values
        local_weights = self._calculate_entropy_weights(local_matrix)
        
        passed['Local_Score'] = self._run_topsis(passed, local_criteria, local_weights)
        
        # Store for first output
        self.calculated_weights = {'Stage 1 (Local)': dict(zip(local_criteria, local_weights))}

        # --- Stage 2: External GTEx Refinement ---
        if ENABLE_EXTERNAL_SPECIFICITY:
            print(f"Stage 2: Refining top {TOP_N_FOR_EXTERNAL} unique genes with GTEx specificity...")
            
            # Select top candidates from local stage
            top_genes = passed.sort_values('Local_Score', ascending=False)['Gene'].unique()[:TOP_N_FOR_EXTERNAL]
            external_data = self.external_manager.fetch_batch_specificity(top_genes)
            
            # Map burden to all rows (gene symbols might repeat across cell types)
            passed['off_target_burden'] = passed['Gene'].map(external_data).fillna(np.nan)
            
            # Recalculate weights including burden for the refined set
            # We only re-calculate for the ones we have burden data for
            refined_subset = passed[passed['off_target_burden'].notna()].copy()
            other_subset = passed[passed['off_target_burden'].isna()].copy()
            
            if not refined_subset.empty:
                full_criteria = local_criteria + ['off_target_burden']
                full_matrix = refined_subset[full_criteria].apply(pd.to_numeric).fillna(0).abs().values
                
                # New "Ideal Weights" discovered after adding external data
                refined_weights = self._calculate_entropy_weights(full_matrix)
                self.calculated_weights['Stage 2 (Refined)'] = dict(zip(full_criteria, refined_weights))
                
                refined_subset['Score'] = self._run_topsis(refined_subset, full_criteria, refined_weights)
                
                # For those not in top N, they keep their local score as final score but penalized
                # Or we can just set them lower.
                other_subset['Score'] = other_subset['Local_Score'] * 0.1 # Significant penalty for not being in validation pool
                
                passed = pd.concat([refined_subset, other_subset])
            else:
                passed['Score'] = passed['Local_Score']
        else:
            passed['Score'] = passed['Local_Score']

        # Final Ranking based on Score
        passed = passed.sort_values(by='Score', ascending=False)
        passed['Rank'] = range(1, len(passed) + 1)
        
        dropped['Rank'] = np.nan
        dropped['Score'] = np.nan
        dropped['Local_Score'] = np.nan
        
        return pd.concat([passed, dropped], ignore_index=True)
