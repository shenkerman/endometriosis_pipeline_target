import pandas as pd
import numpy as np
from .config import CELL_TYPES_INFO, RANKING_WEIGHTS

class PipelineDataProcessor:
    def __init__(self, file_path):
        self.file_path = file_path

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
        """
        Calculates weights using the Entropy Weighting Method.
        1. Normalize matrix: p_ij = x_ij / sum(x_ij)
        2. Calculate entropy: e_j = -k * sum(p_ij * ln(p_ij))
        3. Calculate diversification: d_j = 1 - e_j
        4. Calculate weights: w_j = d_j / sum(d_j)
        """
        # Ensure all values are positive for entropy calculation (using abs)
        # Shift slightly to avoid log(0)
        matrix = matrix + 1e-9
        
        # 1. Normalization
        # Sum of each column
        column_sums = matrix.sum(axis=0)
        p_matrix = matrix / column_sums
        
        # 2. Entropy
        # k = 1 / ln(n_samples)
        n = matrix.shape[0]
        k = 1.0 / np.log(n)
        
        # e = -k * sum(p * ln(p))
        # Handle cases where p=0 by replacing p*ln(p) with 0
        p_ln_p = p_matrix * np.log(p_matrix)
        p_ln_p = np.nan_to_num(p_ln_p)
        entropy = -k * p_ln_p.sum(axis=0)
        
        # 3. Diversification
        diversification = 1.0 - entropy
        
        # 4. Weights
        weights = diversification / diversification.sum()
        
        return weights

    def calculate_ranking(self, df):
        """
        Calculates global ranking using the TOPSIS algorithm with 
        automated Entropy Weighting.
        """
        df = df.copy()
        passed = df[df['Status'] == 'PASS'].copy()
        dropped = df[df['Status'] == 'DROPPED'].copy()
        
        if passed.empty:
            dropped['Rank'] = np.nan
            dropped['Score'] = np.nan
            self.calculated_weights = {}
            return dropped

        # --- Automated Weight Calculation ---
        
        # Prepare criteria matrix
        criteria_cols = ['logFC.eueVSctrl', 'logFC.ecpVSctrl', 'logFC.ecoVSctrl', 'logFC.ecpaVSctrl']
        # Filter for columns that actually exist in the data
        actual_cols = [c for c in criteria_cols if c in passed.columns]
        
        # We transform EuE so that "closer to 0" is "larger value" for entropy calculation 
        # (higher values = more information)
        # However, for pure entropy we use absolute values.
        matrix = passed[actual_cols].apply(pd.to_numeric, errors='coerce').fillna(0).abs().values
        
        # Calculate weights based on data distribution
        weights = self._calculate_entropy_weights(matrix)
        
        # Store for reporting
        self.calculated_weights = dict(zip(actual_cols, weights))

        # --- TOPSIS Implementation ---
        
        # 2. Vector Normalization (TOPSIS specific)
        sums_of_squares = (matrix**2).sum(axis=0)
        norm_matrix = matrix / np.sqrt(sums_of_squares + 1e-9)
        
        # 3. Weighted Normalized Decision Matrix
        weighted_matrix = norm_matrix * weights
        
        # 4. Determine Positive-Ideal and Negative-Ideal Solutions
        v_ideal = []
        v_anti_ideal = []
        
        for i, col in enumerate(actual_cols):
            if col == 'logFC.eueVSctrl':
                # For EuE: Ideal is minimum (0), Anti-Ideal is maximum
                v_ideal.append(weighted_matrix[:, i].min())
                v_anti_ideal.append(weighted_matrix[:, i].max())
            else:
                # For Lesions: Ideal is maximum, Anti-Ideal is minimum (0)
                v_ideal.append(weighted_matrix[:, i].max())
                v_anti_ideal.append(weighted_matrix[:, i].min())
                
        v_ideal = np.array(v_ideal)
        v_anti_ideal = np.array(v_anti_ideal)
        
        # 5. Calculate Geometric Distances
        s_ideal = np.sqrt(((weighted_matrix - v_ideal)**2).sum(axis=1))
        s_anti_ideal = np.sqrt(((weighted_matrix - v_anti_ideal)**2).sum(axis=1))
        
        # 6. Calculate TOPSIS Score
        passed['Score'] = s_anti_ideal / (s_ideal + s_anti_ideal)
        
        # Unified Ranking
        passed = passed.sort_values(by='Score', ascending=False)
        passed['Rank'] = range(1, len(passed) + 1)
        
        dropped['Rank'] = np.nan
        dropped['Score'] = np.nan
        return pd.concat([passed, dropped], ignore_index=True)
