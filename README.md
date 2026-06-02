# Endometriosis Biomarker Pipeline

An autonomous, data-driven pipeline for identifying and ranking biomarker candidates for Endometriosis using multi-criteria decision analysis.

## 🚀 Overview

This project implements a sophisticated bioinformatics pipeline that analyzes single-cell RNA-seq data to discover candidate genes with high diagnostic potential. It balances two often-conflicting biological goals:
1.  **Strong Disease Signal:** Maximum differential expression in ectopic lesions.
2.  **High Specificity:** Minimal expression (logFC near 0) in healthy eutopic tissue.

The pipeline utilizes the **Entropy Weighting Method** for objective feature prioritization and the **TOPSIS algorithm** for robust candidate ranking.

## 📁 Project Structure

```text
├── src/                    # Core logic package
│   ├── config.py           # Configuration and constants
│   ├── processor.py        # Data handling and TOPSIS algorithm
│   └── visualizer.py       # Seaborn-based plotting module
├── run_pipeline.py         # Entry point script
├── data.xlsx               # Primary DE dataset (Tan et al. 2022)
├── results/                # Auto-generated outputs (tables & plots)
└── docs/                   # Detailed documentation
```

## 🛠️ Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/endometriosis-biomarker-pipeline.git
    cd endometriosis-biomarker-pipeline
    ```

2.  **Set up a virtual environment:**
    ```bash
    python3 -m venv venv
    source venv/bin/activate  # On Windows use `venv\Scripts\activate`
    ```

3.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

## 📈 Usage

Simply run the master script to process the data and generate all results:

```bash
python3 run_pipeline.py
```

## 📊 Outputs

The pipeline generates several outputs in the `results/` directory:
*   **Tables:**
    *   `tables/all_pipeline_results.csv`: Unified summary of all genes.
    *   `tables/cell_types/`: Individual results for each subpopulation (dS2, Prv-CCL19, etc.).
*   **Plots:**
    *   `plots/ectopic_vs_eutopic_scatter.png`: Global comparison with Top 10 labels.
    *   `plots/cell_types/`: Individual scatter plots for each cell type, labeling top candidates.
    *   `plots/top_100_candidates_heatmap.png`: High-resolution expression profiles.
*   **Documentation:** Detailed mathematical rationale in `docs/RANKING_LOGIC.md`.

## 🧠 Methodology

For an end-to-end explanation of the filtering (F1, F2, F4) and the ranking algorithms, please see the [Pipeline Summary](docs/PIPELINE_SUMMARY.md) and [Ranking Logic](docs/RANKING_LOGIC.md).

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🔬 Citation

If you use this pipeline in your research, please cite the primary data source:
*Tan et al. (2022). A single-cell atlas of endometriosis. Nature Genetics.*
