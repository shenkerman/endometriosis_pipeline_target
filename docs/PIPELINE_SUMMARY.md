# Endometriosis Biomarker Pipeline: End-to-End Summary

## 1. Project Overview
This pipeline was designed to autonomously identify and rank candidate genes as biomarkers for Endometriosis. It utilizes published single-cell RNA-seq data (Tan et al. 2022) to find genes that are strongly expressed in disease-affected tissues while maintaining high specificity (minimal expression) in healthy eutopic tissue.

The process has been transformed from a manual spreadsheet analysis into a **modular, data-driven software package** that uses advanced multi-criteria decision-making algorithms to select the most promising candidates.

---

## 2. The High-Level Process
The pipeline follows four main stages:
1.  **Data Ingestion:** Parsing pseudo-bulk differential expression (DE) results from `data.xlsx` for 10 distinct cell-type subpopulations.
2.  **Filtering (QC):** Applying rigorous statistical and biological filters to remove noise and inconsistent signals.
3.  **Autonomous Weighting:** Using Information Theory to determine which tissue comparisons provide the most "useful" information for ranking.
4.  **TOPSIS Ranking:** Mathematically scoring each gene based on its geometric distance to a theoretical "Ideal Biomarker."

---

## 3. Detailed Technical Components

### 3.1 Data Source & Scope
*   **Primary Dataset:** Tan et al. 2022.
*   **Tissue Groups:** Eutopic endometrium (EuE), Ectopic peritoneal lesion (EcP), Ectopic ovarian lesion (EcO), and Adjacent peritoneum (EcPA).
*   **Comparison:** All disease groups are compared against a healthy control group (Ctrl).

### 3.2 The Quality & Specificity Filters
Each gene must pass three primary filters to be considered a candidate:
*   **F1 (Noise Floor):** `logCPM ≥ 1`. Prevents inflated Fold Change estimates from low-count genes.
*   **F2 (Significance):** `FDR ≤ 0.01`. Ensures the differential expression is statistically robust.
*   **F4 (Consensus):** The gene must be consistently up-regulated (or down-regulated) in **all** comparisons. This prevents picking genes that are only active in a single tissue type.

### 3.3 Ideal Weight Search (Entropy Weighting Method)
Instead of relying on arbitrary user-defined weights, the pipeline calculates **objective weights** by analyzing the data's distribution:
*   It calculates the **Information Entropy** of each tissue comparison column.
*   Columns with high "Diversification" (where candidates show clear, distinct differences) are automatically assigned higher weights.
*   This ensures the final ranking is purely data-driven and prioritizes the most discriminatory features.

### 3.4 Global Ranking (TOPSIS Algorithm)
The pipeline implements the **TOPSIS** (Technique for Order of Preference by Similarity to Ideal Solution) algorithm.
*   **The Ideal Biomarker ($V^+$):** Defined as having the maximum observed signal in lesions and **exactly zero** signal in the eutopic tissue.
*   **The Anti-Ideal Biomarker ($V^-$):** Defined as having zero signal in lesions and maximum signal in eutopic tissue.
*   **The Score:** Genes are ranked by their "Relative Closeness" to the Ideal. A score of 1.0 would be a perfect biomarker.

---

## 4. Pipeline Outputs & Results
All results are automatically generated in the `results/` directory:

### 4.1 Data Tables (`results/tables/`)
*   **`all_pipeline_results.csv`**: The primary "global heatmap" containing every gene analyzed, their status, and their TOPSIS score.
*   **`candidates_upregulation.csv`**: The top-ranked genes for blood/tissue detection (up-regulated).
*   **`ranking_weights.csv`**: The specific weights determined by the Entropy Weighting Method for that specific run.

### 4.2 Visualizations (`results/plots/`)
*   **Ectopic vs Eutopic Scatter:** Visualizes the trade-off between lesion signal and eutopic specificity, labeling the Top 10 genes.
*   **Status Counts:** Shows the efficiency of the pipeline filters across different cell types.
*   **Top 100 Heatmap:** A high-resolution profile of the 100 most promising biomarkers.

---

## 5. Mathematical Rationale
For a deep dive into the formulas used for Vector Normalization, Euclidean Distance, and Entropy calculation, please refer to the specialized document: **`results/RANKING_LOGIC.md`**.
