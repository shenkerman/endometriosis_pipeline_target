import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from .config import PLOTS_DIR

class PipelineVisualizer:
    def __init__(self, df):
        self.df = df
        self.passed = df[df['Status'] == 'PASS'].copy()
        if not self.passed.empty:
            self.passed['Rank_num'] = pd.to_numeric(self.passed['Rank'], errors='coerce')

    def generate_all_plots(self):
        print("Generating visualizations...")
        sns.set_theme(style="whitegrid")
        os.makedirs(PLOTS_DIR, exist_ok=True)
        
        self.plot_ectopic_vs_eutopic()
        self.plot_status_counts()
        self.plot_top_10_heatmap()
        self.plot_top_100_heatmap()
        self.plot_individual_cell_types()

    def plot_individual_cell_types(self):
        """Generates a separate heatmap for every cell type showing top candidates."""
        if self.passed.empty: return
        
        cell_types_plots_dir = os.path.join(PLOTS_DIR, "cell_types")
        os.makedirs(cell_types_plots_dir, exist_ok=True)
        
        logfc_cols = [c for c in self.passed.columns if isinstance(c, str) and 'logFC' in c]
        
        for cell_type in self.passed['Cell Type'].unique():
            cell_df = self.passed[self.passed['Cell Type'] == cell_type].copy()
            if cell_df.empty: continue
            
            # Take top 20 candidates for this cell type for the heatmap
            top_n_cell = cell_df.sort_values('Rank_num').head(20)
            
            heatmap_data = top_n_cell.set_index('Gene')[logfc_cols].apply(pd.to_numeric).fillna(0)
            
            plt.figure(figsize=(10, 0.5 * len(heatmap_data) + 2))
            sns.heatmap(heatmap_data, annot=True, cmap='RdBu_r', center=0, cbar_kws={'label': 'logFC'})
            
            plt.title(f'Top Candidates logFC Profile: {cell_type}', fontsize=15)
            plt.tight_layout()
            
            # Sanitize filename
            safe_name = cell_type.replace('/', '_').replace(' ', '_')
            plt.savefig(os.path.join(cell_types_plots_dir, f'heatmap_{safe_name}.png'))
            plt.close()

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
        plt.savefig(os.path.join(PLOTS_DIR, 'all_cell_types_scatter.png'))
        plt.close()

    def plot_status_counts(self):
        plt.figure(figsize=(10, 6))
        counts = self.df.groupby(['Cell Type', 'Status']).size().reset_index(name='Count')
        sns.barplot(data=counts, x='Cell Type', y='Count', hue='Status')
        plt.title('Pipeline Status Counts (All Cell Types)', fontsize=16)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(PLOTS_DIR, 'all_cell_types_status_counts.png'))
        plt.close()

    def _generate_heatmap(self, n_top, filename, show_annot=True):
        if self.passed.empty: return
        top_n = self.passed.nsmallest(n_top, 'Rank_num').copy()
        logfc_cols = [c for c in top_n.columns if isinstance(c, str) and 'logFC' in c]
        
        top_n['Label'] = top_n['Gene'] + " (" + top_n['Cell Type'] + ")"
        heatmap_data = top_n.set_index('Label')[logfc_cols].apply(pd.to_numeric).fillna(0)
        
        plt.figure(figsize=(12, 0.5 * min(n_top, 50)))
        sns.heatmap(heatmap_data, annot=show_annot, cmap='RdBu_r', center=0)
        plt.title(f'logFC Profile of Top {n_top} Candidates (All Cell Types)', fontsize=16)
        plt.tight_layout()
        plt.savefig(os.path.join(PLOTS_DIR, filename))
        plt.close()

    def plot_top_10_heatmap(self):
        self._generate_heatmap(10, 'all_cell_types_top_10_heatmap.png', show_annot=True)

    def plot_top_100_heatmap(self):
        self._generate_heatmap(100, 'all_cell_types_top_100_heatmap.png', show_annot=False)
