import os
import pandas as pd
import numpy as np
from src.config import TABLES_DIR, OUTPUT_DIR
from src.processor import PipelineDataProcessor
from src.visualizer import PipelineVisualizer

def main():
    input_file = 'data.xlsx'
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    # 1. Processing
    print("Initializing Modular Pipeline...")
    processor = PipelineDataProcessor(input_file)
    df = processor.load_raw_data()
    df = processor.apply_filters(df)
    df = processor.calculate_ranking(df)

    # 2. Visualizing
    visualizer = PipelineVisualizer(df)
    visualizer.generate_all_plots()

    # 3. Saving Results
    print(f"Saving tables to {TABLES_DIR}...")
    os.makedirs(TABLES_DIR, exist_ok=True)
    
    # Save Automatically Calculated Weights
    weights_df = pd.DataFrame(list(processor.calculated_weights.items()), columns=['Criteria', 'Automated_Weight'])
    weights_df.to_csv(os.path.join(TABLES_DIR, 'all_cell_types_ranking_weights.csv'), index=False)
    
    logfc_cols = sorted([c for c in df.columns if isinstance(c, str) and 'logFC' in c])
    display_cols = ['Rank', 'Gene', 'Cell Type'] + logfc_cols + ['logCPM', 'FDR', 'Score', 'Status', 'Dropped Reason']
    
    # Format Rank for display
    df_display = df.copy()
    df_display['Rank'] = df_display['Rank'].apply(lambda x: f"{int(x)}" if pd.notna(x) else "")
    df_display = df_display[display_cols]
    
    # Save unified results
    df_display.to_csv(os.path.join(TABLES_DIR, 'all_cell_types_results_summary.csv'), index=False)
    
    # Save split results for passed candidates
    passed_mask = df_display['Status'] == 'PASS'
    logfc_ref = logfc_cols[0]
    
    up = df_display[passed_mask & (df[logfc_ref] > 0)].sort_values('Score', ascending=False)
    down = df_display[passed_mask & (df[logfc_ref] < 0)].sort_values('Score', ascending=False)
    
    up.to_csv(os.path.join(TABLES_DIR, 'all_cell_types_candidates_upregulation.csv'), index=False)
    down.to_csv(os.path.join(TABLES_DIR, 'all_cell_types_candidates_downregulation.csv'), index=False)

    # Save Individual Cell Type Tables
    cell_types_tables_dir = os.path.join(TABLES_DIR, "cell_types")
    os.makedirs(cell_types_tables_dir, exist_ok=True)
    
    for cell_type in df_display['Cell Type'].unique():
        cell_df = df_display[df_display['Cell Type'] == cell_type].copy()
        # Sort by Rank (Score)
        cell_df['Rank_int'] = pd.to_numeric(cell_df['Rank'], errors='coerce').fillna(999999)
        cell_df = cell_df.sort_values('Rank_int').drop(columns=['Rank_int'])
        
        safe_name = cell_type.replace('/', '_').replace(' ', '_')
        cell_df.to_csv(os.path.join(cell_types_tables_dir, f'results_{safe_name}.csv'), index=False)

    print("\n" + "="*30)
    print("Modular Pipeline Execution Successful")
    print(f"Results available in: {OUTPUT_DIR}/")
    print("="*30)

if __name__ == "__main__":
    main()
