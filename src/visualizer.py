import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from .config import PLOTS_DIR

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

    def plot_individual_cell_types(self):
        """Generates a separate heatmap for every cell type showing top candidates."""
        if self.passed.empty: return
        
        cell_types_plots_dir = os.path.join(self.output_dir, "cell_types")
        os.makedirs(cell_types_plots_dir, exist_ok=True)
        
        for cell_type in self.passed['Cell Type'].unique():
            cell_df = self.passed[self.passed['Cell Type'] == cell_type].copy()
            if cell_df.empty: continue
            
            # Temporary replace passed with cell_df to use _generate_heatmap
            original_passed = self.passed
            self.passed = cell_df
            
            # Sanitize filename
            safe_name = cell_type.replace('/', '_').replace(' ', '_')
            self._generate_heatmap(20, f'cell_types/heatmap_{safe_name}.png', show_annot=True)
            
            self.passed = original_passed

    def plot_ectopic_vs_eutopic(self):
        if self.passed.empty: return
        plt.figure(figsize=(12, 8))
        ec_cols = ['logFC.ecpVSctrl', 'logFC.ecoVSctrl']
        self.passed['Ectopic_Mean'] = self.passed[ec_cols].mean(axis=1)
        
        sns.scatterplot(data=self.passed, x='logFC.eueVSctrl', y='Ectopic_Mean', 
                        hue='Cell Type', size='logCPM', sizes=(50, 400), alpha=0.7)
        
        top_10 = self.passed.nsmallest(10, 'Rank_num')
        for _, row in top_10.iterrows():
            plt.text(row['logFC.eueVSctrl'], row['Ectopic_Mean'], row['Gene'], 
                     fontsize=9, fontweight='bold', ha='right', va='bottom')
            
        plt.title('Ectopic vs Eutopic DE (Top 10 Candidates Labeled - All Cell Types)', fontsize=16)
        plt.axvline(0, color='grey', linestyle='--', alpha=0.5)
        plt.axhline(0, color='grey', linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'all_cell_types_scatter.png'))
        plt.close()

    def plot_status_counts(self):
        plt.figure(figsize=(10, 6))
        counts = self.df.groupby(['Cell Type', 'Status']).size().reset_index(name='Count')
        sns.barplot(data=counts, x='Cell Type', y='Count', hue='Status')
        plt.title('Pipeline Status Counts (All Cell Types)', fontsize=16)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'all_cell_types_status_counts.png'))
        plt.close()

    def plot_specificity_validation(self):
        """Visualizes the relationship between Ectopic signal and Off-target burden."""
        if self.passed.empty or 'off_target_burden' not in self.passed.columns:
            return

        ec_cols = ['logFC.ecpVSctrl', 'logFC.ecoVSctrl']
        self.passed['Ectopic_Mean'] = self.passed[ec_cols].mean(axis=1)
        
        # Filter for only those that have external data (not 0.0 unless they are truly 0)
        plot_df = self.passed[self.passed['off_target_burden'] >= 0].copy()
        
        fig, ax = plt.subplots(figsize=(12, 8))
        scatter = sns.scatterplot(
            data=plot_df,
            x='logFC.eueVSctrl',
            y='Ectopic_Mean',
            hue='off_target_burden',
            palette='viridis_r', 
            size='Score',
            sizes=(50, 500),
            alpha=0.7,
            edgecolor='black',
            ax=ax
        )
        
        # Label top 10 global candidates
        top_10 = plot_df.nsmallest(10, 'Rank_num')
        for _, row in top_10.iterrows():
            ax.text(row['logFC.eueVSctrl'], row['Ectopic_Mean'], row['Gene'], 
                     fontsize=10, fontweight='bold')
            
        ax.set_title('Integrated Specificity Validation (GTEx + local logFC)', fontsize=16)
        ax.set_xlabel('Eutopic signal (lower is better)', fontsize=12)
        ax.set_ylabel('Ectopic signal (higher is better)', fontsize=12)
        ax.axvline(0, color='grey', linestyle='--', alpha=0.5)
        
        norm = plt.Normalize(plot_df['off_target_burden'].min(), plot_df['off_target_burden'].max())
        sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label='Off-target Burden (GTEx Median TPM)')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'all_cell_types_specificity_validation.png'))
        plt.close()

    def _generate_heatmap(self, n_top, filename, show_annot=True):
        if self.passed.empty: return
        top_n = self.passed.nsmallest(n_top, 'Rank_num').copy()
        
        logfc_cols = [c for c in top_n.columns if isinstance(c, str) and 'logFC' in c]
        has_burden = 'off_target_burden' in top_n.columns and top_n['off_target_burden'].notna().any()
        
        top_n['Label'] = top_n['Gene'] + " (" + top_n['Cell Type'] + ")"
        
        if has_burden:
            # Dual Heatmap Approach
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 0.5 * min(n_top, 50) + 3), 
                                         gridspec_kw={'width_ratios': [len(logfc_cols), 1.2]})
            
            # 1. logFC Heatmap (Diverging)
            logfc_data = top_n.set_index('Label')[logfc_cols].apply(pd.to_numeric).fillna(0)
            sns.heatmap(logfc_data, annot=show_annot, cmap='RdBu_r', center=0, ax=ax1, 
                        cbar_kws={'label': 'logFC (Study)'})
            ax1.set_title('Study Differential Expression', fontsize=14)
            ax1.set_ylabel('')
            
            # 2. GTEx Burden Heatmap (Sequential Spectrum)
            burden_data = top_n.set_index('Label')[['off_target_burden']].apply(pd.to_numeric).fillna(0)
            sns.heatmap(burden_data, annot=show_annot, cmap='magma_r', ax=ax2, 
                        cbar_kws={'label': 'TPM (GTEx Burden)'}, yticklabels=False)
            ax2.set_title('Global Specificity', fontsize=14)
            ax2.set_ylabel('')
            
            plt.suptitle(f'Integrated Biomarker Profile: Top {n_top} Candidates', fontsize=18, y=0.98)
        else:
            # Standard single heatmap
            heatmap_data = top_n.set_index('Label')[logfc_cols].apply(pd.to_numeric).fillna(0)
            plt.figure(figsize=(12, 0.5 * min(n_top, 50) + 2))
            sns.heatmap(heatmap_data, annot=show_annot, cmap='RdBu_r', center=0)
            plt.title(f'logFC Profile of Top {n_top} Candidates', fontsize=16)
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, filename))
        plt.close()

    def plot_top_10_heatmap(self):
        self._generate_heatmap(10, 'all_cell_types_top_10_heatmap.png', show_annot=True)

    def plot_top_100_heatmap(self):
        self._generate_heatmap(100, 'all_cell_types_top_100_heatmap.png', show_annot=False)
