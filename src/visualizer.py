import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from .config import PLOTS_DIR

# Human-readable column name mapping for display in plots
COL_LABELS = {
    'logFC.eueVSctrl':  'EuE',
    'logFC.ecpVSctrl':  'EcP',
    'logFC.ecoVSctrl':  'EcO',
    'logFC.ecpaVSctrl': 'EcPA',
    'specificity_ecp':  '|EcP - EuE|',
    'specificity_eco':  '|EcO - EuE|',
    'gtex_burden':      'GTEx off-target',
    'off_target_burden':'GTEx off-target',
    'cellxgene_burden': 'CellxGene off-target',
}


class PipelineVisualizer:
    def __init__(self, df, output_dir=None):
        self.df = df
        self.output_dir = output_dir if output_dir else PLOTS_DIR
        self.passed = df[df['Status'] == 'PASS'].copy()
        if not self.passed.empty:
            self.passed['Rank_num'] = pd.to_numeric(self.passed['Rank'], errors='coerce')

    def generate_all_plots(self):
        print("Generating visualizations...")
        sns.set_theme(style="whitegrid")
        os.makedirs(self.output_dir, exist_ok=True)

        self.plot_ectopic_vs_eutopic()
        self.plot_status_counts()
        self.plot_top_10_heatmap()
        self.plot_top_100_heatmap()
        self.plot_individual_cell_types()
        self.plot_specificity_validation()

    # ------------------------------------------------------------------ #
    #  Individual cell type heatmaps
    # ------------------------------------------------------------------ #
    def plot_individual_cell_types(self):
        if self.passed.empty:
            return
        cell_types_plots_dir = os.path.join(self.output_dir, "cell_types")
        os.makedirs(cell_types_plots_dir, exist_ok=True)

        for cell_type in self.passed['Cell Type'].unique():
            cell_df = self.passed[self.passed['Cell Type'] == cell_type].copy()
            if cell_df.empty:
                continue
            original_passed = self.passed
            self.passed = cell_df
            safe_name = cell_type.replace('/', '_').replace(' ', '_')
            self._generate_heatmap(20, f'cell_types/heatmap_{safe_name}.png', show_annot=True)
            self.passed = original_passed

    # ------------------------------------------------------------------ #
    #  Scatter: Ectopic vs Eutopic (overview)
    # ------------------------------------------------------------------ #
    def plot_ectopic_vs_eutopic(self):
        if self.passed.empty:
            return
        plt.figure(figsize=(12, 8))
        ec_cols = [c for c in ['logFC.ecpVSctrl', 'logFC.ecoVSctrl'] if c in self.passed.columns]
        if not ec_cols:
            return
        self.passed['Ectopic_Mean'] = self.passed[ec_cols].mean(axis=1)

        sns.scatterplot(
            data=self.passed,
            x='logFC.eueVSctrl',
            y='Ectopic_Mean',
            hue='Cell Type',
            size='logCPM',
            sizes=(50, 400),
            alpha=0.7
        )

        top_10 = self.passed.nsmallest(10, 'Rank_num')
        for _, row in top_10.iterrows():
            plt.text(row['logFC.eueVSctrl'], row['Ectopic_Mean'], row['Gene'],
                     fontsize=9, fontweight='bold', ha='right', va='bottom')

        plt.title('Ectopic vs Eutopic Expression — All Cell Types\n(Top 10 candidates labeled)',
                  fontsize=15)
        plt.xlabel('EuE logFC (eutopic endometrium vs healthy control)', fontsize=12)
        plt.ylabel('Mean Lesion logFC — mean(EcP, EcO) vs healthy control', fontsize=12)
        plt.axvline(0, color='grey', linestyle='--', alpha=0.5)
        plt.axhline(0, color='grey', linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'all_cell_types_scatter.png'))
        plt.close()

    # ------------------------------------------------------------------ #
    #  Bar: filter status counts
    # ------------------------------------------------------------------ #
    def plot_status_counts(self):
        plt.figure(figsize=(10, 6))
        counts = self.df.groupby(['Cell Type', 'Status']).size().reset_index(name='Count')
        sns.barplot(data=counts, x='Cell Type', y='Count', hue='Status')
        plt.title('Pipeline Filter Status — All Cell Types', fontsize=15)
        plt.xlabel('Cell Type', fontsize=12)
        plt.ylabel('Number of Genes', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'all_cell_types_status_counts.png'))
        plt.close()

    # ------------------------------------------------------------------ #
    #  Scatter: Lesion Specificity vs Off-Target Burden (V2)
    # ------------------------------------------------------------------ #
    def plot_specificity_validation(self):
        """
        V2: X = |EcP - EuE| (peritoneal lesion specificity)
            Y = |EcO - EuE| (ovarian lesion specificity)
            Color = GTEx off-target burden (log2 ratio vs uterus)
            Size = TOPSIS Score
        Top-right quadrant + green = ideal candidate.
        """
        if self.passed.empty:
            return

        # Use V2 specificity columns if available, fall back to V1 raw logFC
        has_specificity = ('specificity_ecp' in self.passed.columns and
                           'specificity_eco' in self.passed.columns)

        if has_specificity:
            plot_df = self.passed.copy()
            x_col   = 'specificity_ecp'
            y_col   = 'specificity_eco'
            x_label = '|EcP − EuE|  (peritoneal lesion specificity, log₂FC)'
            y_label = '|EcO − EuE|  (ovarian lesion specificity, log₂FC)'
            threshold = 1.0   # 2-fold in log2 space
        else:
            # Fall back to V1 style (no specificity columns computed)
            ec_cols = [c for c in ['logFC.ecpVSctrl', 'logFC.ecoVSctrl'] if c in self.passed.columns]
            if not ec_cols:
                return
            self.passed['Ectopic_Mean'] = self.passed[ec_cols].mean(axis=1)
            plot_df = self.passed.copy()
            x_col   = 'logFC.eueVSctrl'
            y_col   = 'Ectopic_Mean'
            x_label = 'EuE logFC (eutopic endometrium vs control)'
            y_label = 'Mean Lesion logFC — mean(EcP, EcO)'
            threshold = None

        # Determine which off-target column to use for color
        burden_col = None
        burden_label = ''
        for col, label in [
            ('gtex_burden',      'Off-target burden\nmax log₂(GTEx / uterus TPM)'),
            ('off_target_burden','Off-target burden\nmax log₂(GTEx / uterus TPM)'),
        ]:
            if col in plot_df.columns and plot_df[col].notna().any():
                burden_col  = col
                burden_label = label
                break

        # Only keep rows with off-target data if available, otherwise show all
        if burden_col:
            plot_df = plot_df[plot_df[burden_col].notna()].copy()
        if plot_df.empty:
            return

        fig, ax = plt.subplots(figsize=(12, 9))

        scatter_kwargs = dict(
            data=plot_df,
            x=x_col,
            y=y_col,
            size='Score',
            sizes=(60, 500),
            alpha=0.75,
            edgecolor='black',
            linewidth=0.4,
            ax=ax
        )

        if burden_col:
            scatter_kwargs['hue']     = burden_col
            scatter_kwargs['palette'] = 'RdYlGn_r'   # Green = low off-target (good)
        else:
            scatter_kwargs['hue'] = 'Cell Type'

        sns.scatterplot(**scatter_kwargs)

        # Add colorbar for off-target burden
        if burden_col:
            norm = plt.Normalize(plot_df[burden_col].min(), plot_df[burden_col].max())
            sm   = plt.cm.ScalarMappable(cmap='RdYlGn_r', norm=norm)
            sm.set_array([])
            fig.colorbar(sm, ax=ax, label=burden_label, pad=0.02)
            ax.get_legend().remove()

        # Label top candidates
        top_n = min(10, len(plot_df))
        top_genes = plot_df.nsmallest(top_n, 'Rank_num')
        for _, row in top_genes.iterrows():
            ax.text(row[x_col], row[y_col], row['Gene'],
                    fontsize=9, fontweight='bold', ha='left', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.6, ec='none'))

        # Reference lines
        ax.axvline(0, color='grey', linestyle='--', alpha=0.4, linewidth=0.8)
        ax.axhline(0, color='grey', linestyle='--', alpha=0.4, linewidth=0.8)
        if threshold is not None:
            ax.axvline(threshold, color='steelblue', linestyle=':', alpha=0.5,
                       linewidth=1, label=f'2-fold threshold (log₂FC = {threshold})')
            ax.axhline(threshold, color='steelblue', linestyle=':', alpha=0.5, linewidth=1)

        ax.set_title(
            'Lesion Specificity vs Off-Target Burden — Top Candidates\n'
            'Top-right + green = high specificity in both lesion types, low off-target risk',
            fontsize=14
        )
        ax.set_xlabel(x_label, fontsize=12)
        ax.set_ylabel(y_label, fontsize=12)

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'all_cell_types_specificity_validation.png'),
                    dpi=150, bbox_inches='tight')
        plt.close()

    # ------------------------------------------------------------------ #
    #  Heatmap (shared)
    # ------------------------------------------------------------------ #
    def _generate_heatmap(self, n_top, filename, show_annot=True):
        if self.passed.empty:
            return
        top_n = self.passed.nsmallest(n_top, 'Rank_num').copy()

        # Use specificity columns if available, otherwise raw logFC
        spec_cols  = [c for c in ['specificity_ecp', 'specificity_eco'] if c in top_n.columns]
        logfc_cols = [c for c in top_n.columns if isinstance(c, str) and 'logFC' in c]

        display_cols = spec_cols if spec_cols else logfc_cols
        has_burden   = 'off_target_burden' in top_n.columns and top_n['off_target_burden'].notna().any()

        top_n['Label'] = top_n['Gene'] + '  (' + top_n['Cell Type'] + ')'

        # Rename columns for display
        rename_map = {c: COL_LABELS.get(c, c) for c in display_cols}

        if has_burden:
            fig, (ax1, ax2) = plt.subplots(
                1, 2, figsize=(16, 0.5 * min(n_top, 50) + 3),
                gridspec_kw={'width_ratios': [len(display_cols), 1.2]}
            )
            data1 = top_n.set_index('Label')[display_cols].rename(columns=rename_map).apply(pd.to_numeric).fillna(0)
            sns.heatmap(data1, annot=show_annot, fmt='.2f', cmap='YlOrRd', ax=ax1,
                        cbar_kws={'label': 'Specificity distance (|logFC|)'})
            ax1.set_title('Lesion Specificity (log₂FC)', fontsize=13)
            ax1.set_ylabel('')

            data2 = top_n.set_index('Label')[['off_target_burden']].rename(
                columns={'off_target_burden': 'GTEx off-target\n(log₂ ratio vs uterus)'}
            ).apply(pd.to_numeric).fillna(0)
            sns.heatmap(data2, annot=show_annot, fmt='.2f', cmap='magma_r', ax=ax2,
                        cbar_kws={'label': 'log₂ off-target ratio'}, yticklabels=False)
            ax2.set_title('Off-Target Burden\n(GTEx, uterus-normalized)', fontsize=13)
            ax2.set_ylabel('')

            plt.suptitle(
                f'Integrated Biomarker Profile — Top {n_top} Candidates\n'
                'Source: Tan et al. 2022 (Supp. Table 5) | Off-target: GTEx v8',
                fontsize=15, y=1.01
            )
        else:
            heatmap_data = top_n.set_index('Label')[display_cols].rename(columns=rename_map).apply(pd.to_numeric).fillna(0)
            plt.figure(figsize=(10, 0.5 * min(n_top, 50) + 2))
            sns.heatmap(heatmap_data, annot=show_annot, fmt='.2f', cmap='YlOrRd',
                        cbar_kws={'label': 'Specificity distance (|logFC|)'})
            plt.title(
                f'Lesion Specificity Profile — Top {n_top} Candidates\n'
                'Source: Tan et al. 2022 (Supp. Table 5)',
                fontsize=14
            )

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, filename), dpi=150, bbox_inches='tight')
        plt.close()

    def plot_top_10_heatmap(self):
        self._generate_heatmap(10, 'all_cell_types_top_10_heatmap.png', show_annot=True)

    def plot_top_100_heatmap(self):
        self._generate_heatmap(100, 'all_cell_types_top_100_heatmap.png', show_annot=False)
