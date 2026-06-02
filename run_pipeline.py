import os
import datetime
import pandas as pd
import numpy as np
from src.config import TABLES_DIR, OUTPUT_DIR, RUNS_DIR, ENABLE_EXTERNAL_SPECIFICITY
from src.processor import PipelineDataProcessor
from src.visualizer import PipelineVisualizer

def main():
    input_file = 'data.xlsx'
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    # 0. Setup Run Directory
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_id = f"RUN_{timestamp}"
    current_run_dir = os.path.join(RUNS_DIR, run_id)
    run_tables_dir = os.path.join(current_run_dir, "tables")
    run_plots_dir = os.path.join(current_run_dir, "plots")
    
    os.makedirs(run_tables_dir, exist_ok=True)
    os.makedirs(run_plots_dir, exist_ok=True)

    # 1. Processing
    print(f"--- Starting Pipeline Run: {run_id} ---")
    print(f"External Data Enabled: {ENABLE_EXTERNAL_SPECIFICITY}")
    
    processor = PipelineDataProcessor(input_file)
    df = processor.load_raw_data()
    df = processor.apply_filters(df)
    df = processor.calculate_ranking(df)

    # 2. Visualizing
    visualizer = PipelineVisualizer(df, output_dir=run_plots_dir)
    visualizer.generate_all_plots()

    # 3. Saving Results
    print(f"Saving results to {current_run_dir}...")
    
    # Save Automatically Calculated Weights
    weights_df = pd.DataFrame(list(processor.calculated_weights.items()), columns=['Criteria', 'Automated_Weight'])
    weights_df.to_csv(os.path.join(run_tables_dir, 'all_cell_types_ranking_weights.csv'), index=False)
    
    # Save Dedicated External Specificity Validation Table
    if ENABLE_EXTERNAL_SPECIFICITY:
        validation_df = df[df['off_target_burden'].notna()].copy()
        validation_cols = ['Rank', 'Gene', 'Cell Type', 'Score', 'off_target_burden', 'logFC.eueVSctrl']
        validation_df = validation_df[validation_cols].sort_values('Score', ascending=False)
        validation_df.to_csv(os.path.join(run_tables_dir, 'all_cell_types_external_specificity_validation.csv'), index=False)
    
    logfc_cols = [c for c in df.columns if isinstance(c, str) and 'logFC' in c]
    if 'off_target_burden' in df.columns:
        logfc_cols.append('off_target_burden')
    
    display_cols = ['Rank', 'Gene', 'Cell Type'] + sorted(logfc_cols) + ['logCPM', 'FDR', 'Score', 'Status', 'Dropped Reason']
    
    # Format Rank for display
    df_display = df.copy()
    df_display['Rank'] = df_display['Rank'].apply(lambda x: f"{int(x)}" if pd.notna(x) else "")
    
    # Ensure all display columns exist
    for col in display_cols:
        if col not in df_display.columns:
            df_display[col] = np.nan
            
    df_display = df_display[display_cols]
    
    # Save unified results
    df_display.to_csv(os.path.join(run_tables_dir, 'all_cell_types_results_summary.csv'), index=False)
    
    # Save split results for passed candidates
    passed_mask = df_display['Status'] == 'PASS'
    logfc_ref = [c for c in logfc_cols if 'logFC' in c][0]
    
    up = df_display[passed_mask & (df[logfc_ref] > 0)].sort_values('Score', ascending=False)
    down = df_display[passed_mask & (df[logfc_ref] < 0)].sort_values('Score', ascending=False)
    
    up.to_csv(os.path.join(run_tables_dir, 'all_cell_types_candidates_upregulation.csv'), index=False)
    down.to_csv(os.path.join(run_tables_dir, 'all_cell_types_candidates_downregulation.csv'), index=False)

    # Save Individual Cell Type Tables
    cell_types_tables_dir = os.path.join(run_tables_dir, "cell_types")
    os.makedirs(cell_types_tables_dir, exist_ok=True)
    
    for cell_type in df_display['Cell Type'].unique():
        cell_df = df_display[df_display['Cell Type'] == cell_type].copy()
        cell_df['Rank_int'] = pd.to_numeric(cell_df['Rank'], errors='coerce').fillna(999999)
        cell_df = cell_df.sort_values('Rank_int').drop(columns=['Rank_int'])
        
        safe_name = cell_type.replace('/', '_').replace(' ', '_')
        cell_df.to_csv(os.path.join(cell_types_tables_dir, f'results_{safe_name}.csv'), index=False)

    # 4. Update Global History Log
    history_file = os.path.join(OUTPUT_DIR, 'run_history.csv')
    history_entry = {
        'Run_ID': run_id,
        'Timestamp': timestamp,
        'External_Data': ENABLE_EXTERNAL_SPECIFICITY,
        'Path': current_run_dir,
        'Weights': str(processor.calculated_weights)
    }
    history_df = pd.DataFrame([history_entry])
    if os.path.exists(history_file):
        history_df = pd.concat([pd.read_csv(history_file), history_df], ignore_index=True)
    history_df.to_csv(history_file, index=False)
    
    # 5. Copy docs for completeness
    run_docs_dir = os.path.join(current_run_dir, "docs")
    os.makedirs(run_docs_dir, exist_ok=True)
    import shutil
    for doc_file in os.listdir("docs"):
        if doc_file.endswith(".md"):
            shutil.copy(os.path.join("docs", doc_file), run_docs_dir)

    print("\n" + "="*30)
    print("Pipeline Execution Successful")
    print(f"Run Directory: {current_run_dir}")
    print("="*30)

if __name__ == "__main__":
    main()
