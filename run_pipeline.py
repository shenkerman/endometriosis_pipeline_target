import os
import argparse
import datetime
import pandas as pd
import numpy as np
import shutil
from src.config import OUTPUT_DIR, RUNS_DIR, ENABLE_EXTERNAL_SPECIFICITY
from src.processor import PipelineDataProcessor
from src.visualizer import PipelineVisualizer


def parse_args():
    parser = argparse.ArgumentParser(
        description='Endometriosis Biomarker Pipeline V2',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 run_pipeline.py                          # all cell types
  python3 run_pipeline.py --cell-type dS2          # dS2 only
  python3 run_pipeline.py --cell-type dS2 dS1      # dS2 and dS1
  python3 run_pipeline.py --cell-type "Prv-CCL19"  # names with hyphens need quotes
        """
    )
    parser.add_argument(
        '--cell-type',
        nargs='+',
        metavar='CELL_TYPE',
        default=None,
        help='One or more cell types to include. If omitted, all cell types are used.'
    )
    return parser.parse_args()


def main():
    args = parse_args()

    input_file = 'data.xlsx'
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    # 0. Setup run directory
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_id = f"RUN_{timestamp}"
    if args.cell_type:
        safe_ct = '_'.join(ct.replace('/', '_').replace(' ', '_') for ct in args.cell_type)
        run_id += f"_CT_{safe_ct}"

    current_run_dir = os.path.join(RUNS_DIR, run_id)
    run_tables_dir  = os.path.join(current_run_dir, "tables")
    run_plots_dir   = os.path.join(current_run_dir, "plots")
    os.makedirs(run_tables_dir, exist_ok=True)
    os.makedirs(run_plots_dir, exist_ok=True)

    # 1. Load and filter
    print(f"\n--- Pipeline V2 | Run: {run_id} ---")
    print(f"External data: {ENABLE_EXTERNAL_SPECIFICITY}")
    print(f"Cell type filter: {args.cell_type if args.cell_type else 'None (all cell types)'}")

    processor = PipelineDataProcessor(input_file)
    df = processor.load_raw_data()

    if args.cell_type:
        df = df[df['Cell Type'].isin(args.cell_type)].copy()
        if df.empty:
            all_types = processor.load_raw_data()['Cell Type'].unique().tolist()
            print(f"Error: No data found for: {args.cell_type}")
            print(f"Available cell types: {all_types}")
            return

    df = processor.apply_filters(df)
    df = processor.calculate_ranking(df)

    # 2. Visualize
    visualizer = PipelineVisualizer(df, output_dir=run_plots_dir)
    visualizer.generate_all_plots()

    # 3. Build display columns
    logfc_cols = sorted([c for c in df.columns if isinstance(c, str) and 'logFC' in c])
    v2_cols    = [c for c in ['specificity_ecp', 'specificity_eco',
                               'gtex_burden', 'cellxgene_burden', 'off_target_agree',
                               'CPM', 'logCPM', 'CPM_percentile', 'EuE_near_zero']
                  if c in df.columns]

    display_cols = (
        ['Rank', 'Gene', 'Cell Type']
        + logfc_cols
        + v2_cols
        + [c for c in ['FDR', 'Score', 'Local_Score', 'Status', 'Dropped Reason'] if c in df.columns]
    )

    df_display = df.copy()
    df_display['Rank'] = df_display['Rank'].apply(lambda x: f"{int(x)}" if pd.notna(x) else "")
    for col in display_cols:
        if col not in df_display.columns:
            df_display[col] = np.nan
    df_display = df_display[display_cols]

    # 4. Save results
    print(f"\nSaving results to {current_run_dir}...")

    # Weights table
    weights_rows = [
        {'Stage': stage, 'Criteria': crit, 'Weight': round(val, 4)}
        for stage, w in processor.calculated_weights.items()
        for crit, val in w.items()
    ]
    pd.DataFrame(weights_rows).to_csv(
        os.path.join(run_tables_dir, 'all_cell_types_ranking_weights.csv'), index=False)

    # Full results
    df_display.to_csv(
        os.path.join(run_tables_dir, 'all_cell_types_results_summary.csv'), index=False)

    # External specificity validation table
    if ENABLE_EXTERNAL_SPECIFICITY and 'gtex_burden' in df.columns:
        val_cols = [c for c in ['Rank', 'Gene', 'Cell Type', 'Score', 'Local_Score',
                                 'specificity_ecp', 'specificity_eco',
                                 'gtex_burden', 'cellxgene_burden', 'off_target_agree',
                                 'logFC.eueVSctrl', 'CPM_percentile', 'EuE_near_zero']
                    if c in df_display.columns]
        df_display[df_display['Status'] == 'PASS'][val_cols].to_csv(
            os.path.join(run_tables_dir, 'all_cell_types_external_specificity_validation.csv'),
            index=False)

    # UP / DOWN split
    # V2: direction determined by mean(EcP, EcO), not just EcO (V1 bug fix)
    passed_mask = df['Status'] == 'PASS'
    ecp_vals    = pd.to_numeric(df.get('logFC.ecpVSctrl', pd.Series(np.nan, index=df.index)), errors='coerce')
    eco_vals    = pd.to_numeric(df.get('logFC.ecoVSctrl', pd.Series(np.nan, index=df.index)), errors='coerce')
    mean_lesion = (ecp_vals + eco_vals) / 2

    # Genes with mean_lesion > 0: upregulated in lesions vs control
    # Genes with mean_lesion <= 0: downregulated or neutral (included in DOWN to avoid silent drops)
    up   = df_display[passed_mask & (mean_lesion > 0)].copy()
    down = df_display[passed_mask & (mean_lesion <= 0)].copy()
    up.to_csv(os.path.join(run_tables_dir,   'all_cell_types_candidates_upregulation.csv'),   index=False)
    down.to_csv(os.path.join(run_tables_dir, 'all_cell_types_candidates_downregulation.csv'), index=False)

    # Per-cell-type tables
    cell_types_dir = os.path.join(run_tables_dir, "cell_types")
    os.makedirs(cell_types_dir, exist_ok=True)
    for ct in df_display['Cell Type'].unique():
        ct_df = df_display[df_display['Cell Type'] == ct].copy()
        ct_df['_rank_int'] = pd.to_numeric(ct_df['Rank'], errors='coerce').fillna(999999)
        ct_df = ct_df.sort_values('_rank_int').drop(columns=['_rank_int'])
        safe  = ct.replace('/', '_').replace(' ', '_')
        ct_df.to_csv(os.path.join(cell_types_dir, f'results_{safe}.csv'), index=False)

    # 5. Run history log
    history_file = os.path.join(OUTPUT_DIR, 'run_history.csv')
    entry = pd.DataFrame([{
        'Run_ID':        run_id,
        'Timestamp':     timestamp,
        'Cell_Type':     str(args.cell_type) if args.cell_type else 'All',
        'External_Data': ENABLE_EXTERNAL_SPECIFICITY,
        'Path':          current_run_dir,
        'Weights':       str(processor.calculated_weights)
    }])
    if os.path.exists(history_file):
        entry = pd.concat([pd.read_csv(history_file), entry], ignore_index=True)
    entry.to_csv(history_file, index=False)

    # 6. Copy docs
    run_docs_dir = os.path.join(current_run_dir, "docs")
    os.makedirs(run_docs_dir, exist_ok=True)
    if os.path.exists("docs"):
        for doc_file in os.listdir("docs"):
            if doc_file.endswith(".md"):
                shutil.copy(os.path.join("docs", doc_file), run_docs_dir)

    print("\n" + "=" * 30)
    print("Pipeline V2 - Execution Successful")
    print(f"Run directory: {current_run_dir}")
    print("=" * 30)


if __name__ == "__main__":
    main()
