# Endometriosis Biomarker Pipeline — V2

An autonomous, data-driven pipeline for identifying and ranking candidate biomarker genes for endometriosis using multi-criteria decision analysis (TOPSIS).

V2 fixes several biological and algorithmic issues identified during manual validation against Tan et al. 2022 (Supplementary Table 5). See [Changes from V1](#changes-from-v1) for details.

---

## Overview

The pipeline analyzes single-cell RNA-seq differential expression data to discover candidate genes with high diagnostic or therapeutic potential. It balances two biological goals:

1. **Lesion specificity:** strong differential expression in ectopic lesions relative to the patient's own eutopic endometrium
2. **Off-target safety:** minimal expression in healthy non-reproductive tissues

### Data source

Tan, Yuliana et al. "Single-cell analysis of endometriosis reveals a coordinated transcriptional programme driving immunotolerance and angiogenesis across eutopic and ectopic tissues." *Nature Cell Biology* 24.8 (2022): 1306–1318. [doi:10.1038/s41556-022-00961-5](https://doi.org/10.1038/s41556-022-00961-5)

Supplementary Table 5: pseudo-bulk differential expression (edgeR glmQLFTest) per cell type across four tissue comparisons: EuE (eutopic endometrium), EcP (ectopic peritoneal lesion), EcO (ectopic ovarian lesion), EcPA (ectopic peritoneal adjacent).

---

## Project Structure

```
├── src/
│   ├── config.py          # All configuration and constants
│   ├── processor.py       # Filtering, TOPSIS ranking, derived columns
│   ├── external_data.py   # GTEx API + CellxGene off-target data
│   └── visualizer.py      # Seaborn-based plotting
├── run_pipeline.py        # Entry point with CLI arguments
├── data.xlsx              # Primary DE dataset (Tan et al. 2022, Supp. Table 5)
├── results/               # Auto-generated outputs (tables and plots)
└── docs/                  # Detailed documentation
```

---

## Installation

```bash
git clone https://github.com/shenkerman/endometriosis_pipeline_target.git
cd endometriosis_pipeline_target
```

**Option A - with CellxGene off-target support (Python 3.12 required):**

```bash
conda create -n endo_pipeline python=3.12 -y
conda activate endo_pipeline
pip install cellxgene-census matplotlib seaborn pandas numpy requests scipy openpyxl
```

**Option B - GTEx-only mode (any Python 3.8+):**

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

> `cellxgene-census` requires Python < 3.13. If you are on Python 3.13+, use Option A with conda.

---

## Usage

Activate the environment first (if using conda):
```bash
conda activate endo_pipeline
```

Run on all cell types (default):
```bash
python3 run_pipeline.py
```

Run on a specific cell type:
```bash
python3 run_pipeline.py --cell-type dS2
```

Run on multiple cell types:
```bash
python3 run_pipeline.py --cell-type dS2 dS1
python3 run_pipeline.py --cell-type "Prv-CCL19"   # use quotes for names with hyphens
```

If `--cell-type` is not provided, all cell types are included with equal weight.

---

## Outputs

All outputs are saved to a timestamped directory under `results/runs/RUN_<timestamp>/`.

**Tables:**
- `all_cell_types_results_summary.csv` — full ranked results with all columns
- `all_cell_types_candidates_upregulation.csv` — passed genes with mean(EcP, EcO) > 0
- `all_cell_types_candidates_downregulation.csv` — passed genes with mean(EcP, EcO) < 0
- `all_cell_types_ranking_weights.csv` — entropy weights used in each TOPSIS stage
- `all_cell_types_external_specificity_validation.csv` — off-target data for validated genes
- `tables/cell_types/results_<cell_type>.csv` — per-cell-type breakdown

**Plots:**
- `plots/ectopic_vs_eutopic_scatter.png` — global scatter with top 10 labels
- `plots/cell_types/` — per-cell-type scatter plots
- `plots/top_100_candidates_heatmap.png` — expression heatmap

**Plot legend (for all plots including off-target axes):**

*Analysis:* Candidate genes ranked by lesion specificity and off-target expression using TOPSIS multi-criteria scoring. Source: Tan et al. 2022 (see above).

*Tissue types:* EuE — eutopic endometrium (patient-matched control); EcP — ectopic peritoneal lesion; EcO — ectopic ovarian lesion. Specificity = |EcP − EuE| or |EcO − EuE| (log2FC).

*Cell types:* dS2 — decidual stromal type 2 (dominant stromal population in ectopic lesions); dS1 — decidual stromal type 1; Prv-CCL19 — perivascular CCL19+ cells; EC-tip — tip endothelial cells; EC-aPCV — arterial post-capillary venule endothelial cells; EC-PCV — post-capillary venule endothelial cells; Mɸ1-LYVE1 — LYVE1+ macrophages; Mɸ4-infiltrated — infiltrated macrophages; cDC2 — conventional dendritic cells type 2; Treg — regulatory T cells; B cell — B lymphocytes.

*Off-target color scale (all plots):* Color = max log2((GTEx TPM + 1) / (uterus TPM + 1)) across non-reproductive tissues. Green = low off-target risk; red = broadly expressed in healthy tissue.

> **Note on bulk vs single-cell comparability:** Off-target axes use GTEx bulk RNA-seq normalized to GTEx "Uterus" tissue (mixed bulk: myometrium + endometrium). This is not equivalent to EuE in Tan et al. 2022, which is scRNA-seq from eutopic endometrial cells. Off-target values are directionally valid but not directly interchangeable with logFC from Supplementary Table 5.

---

## Methodology

### Filters

| Filter | Criterion | Rationale |
|--------|-----------|-----------|
| F1 | logCPM >= 1 | Removes genes below reliable expression threshold (Chen et al. 2016; Squair et al. 2021) |
| F2 | FDR <= 0.01 | Controls false discovery rate (edgeR guide; OSCA) |

F4 ("inconsistent direction") was removed in V2. A gene downregulated in EuE but upregulated in lesions is exactly the desired biomarker profile and should not be penalized.

### Ranking — TOPSIS with Lesion-Specificity Distances

**Stage 1 (local):** TOPSIS on two derived criteria:

```
specificity_EcP = |EcP logFC - EuE logFC|
specificity_EcO = |EcO logFC - EuE logFC|
```

These replace the raw logFC columns used in V1. The distance formulation captures both signal strength and specificity in one number — a gene with EcP = +3 and EuE = 0 and a gene with EcP = −3 and EuE = 0 score equally (distance = 3), correctly treating both directions as valid. EcPA is excluded.

Weights are computed via entropy weighting, applied to the specificity distances (not raw logFC values).

**Stage 2 (external refinement):** For the top 200 genes by local score, two off-target metrics are added as minimize criteria:

| Metric | Source | What it measures |
|--------|--------|-----------------|
| `gtex_burden_vs_uterus` | GTEx API (v8) | max log2((TPM_off-target + 1) / (TPM_uterus + 1)) across non-reproductive tissues |
| `cellxgene_burden` | CellxGene Census API | max % cells expressing the gene across non-reproductive healthy tissues (disease = normal, is_primary_data = True) |

GTEx metric uses uterus as reference so values are normalized per gene's baseline reproductive expression. Negative ratios (gene is less expressed off-target than in uterus) are clipped to 0.

Genes outside the top 200 retain their Stage 1 local score with no penalty (V1 applied a 0.1x multiplier, which has been removed).

`off_target_agree` column: True if GTEx and CellxGene agree on risk level (both above or both below their respective dataset medians). False = disagreement, flag for manual review.

### UP / DOWN Classification

Output tables are split into upregulated and downregulated candidates based on:

```
mean(EcP logFC, EcO logFC) > 0   ->  UP
mean(EcP logFC, EcO logFC) <= 0  ->  DOWN  (includes genes at exactly 0)
```

V1 used only the EcO column (alphabetically first), which incorrectly classified genes with strong EcP signal but weak EcO signal.

### Additional Output Columns

| Column | Description |
|--------|-------------|
| `specificity_ecp` | \|EcP - EuE\| logFC |
| `specificity_eco` | \|EcO - EuE\| logFC |
| `gtex_burden_vs_uterus` | max log₂(TPM in off-target tissue / TPM in uterus) — clipped at 0; uterus-normalized GTEx off-target burden |
| `CPM` | Back-transformed expression: 2^logCPM |
| `CPM_percentile` | Rank percentile of CPM within each cell type (top 15% sorted first) |
| `off_target_agree` | True if GTEx and CellxGene agree on risk level |

Results are sorted with top-15%-by-CPM genes first within each score tier (high-confidence expression signal prioritized).

---

## CellxGene Data

CellxGene off-target integration uses the [CellxGene Census Python API](https://chanzuckerberg.github.io/cellxgene-census/) to fetch % cells expressing each gene across healthy non-reproductive tissues (`disease == 'normal'`, `is_primary_data == True`). Data is fetched in a single batch for all top-N candidates and cached locally in `.api_cache/cellxgene_cache.pkl`.

**Requirement:** `cellxgene-census` requires **Python < 3.13**. Install via conda:

```bash
conda create -n endo_pipeline python=3.12 -y
conda activate endo_pipeline
pip install cellxgene-census
```

If `cellxgene-census` is not installed, CellxGene integration is skipped automatically and only GTEx is used. The first Census fetch may take several minutes; subsequent runs use the local cache.

---

## Reference Specifications

| Method | Reference |
|--------|-----------|
| Pseudo-bulk DE (edgeR glmQLFTest) | Squair et al. 2021, *Nature Communications* — [doi:10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2) |
| logCPM >= 1 noise floor | Chen et al. 2016, *F1000Research* — [doi:10.12688/f1000research.9005.2](https://doi.org/10.12688/f1000research.9005.2) |
| Pseudo-bulk workflow | OSCA — [bioconductor.org/books/release/OSCA](https://bioconductor.org/books/release/OSCA) |
| edgeR glmQLFTest | edgeR User's Guide — [bioconductor.org/packages/edgeR](https://bioconductor.org/packages/edgeR) |
| TOPSIS | Hwang & Yoon 1981, *Multiple Attribute Decision Making* |
| Entropy weighting | Shannon 1948; applied per column to reward discriminating criteria |

---

## Changes from V1

| # | Change | Impact |
|---|--------|--------|
| 1 | Cell type runtime parameter (`--cell-type`) | Focus analysis on biologically relevant populations without editing code |
| 2 | CellxGene Census API off-target integration (replaces manual CSV export) | Adds breadth metric alongside GTEx intensity metric; fully automated via Python API |
| 3 | GTEx metric: uterus-normalized log2 ratio (replaces raw TPM sum) | Comparable across genes; penalizes spikes over uterus baseline, not cumulative breadth |
| 4 | `off_target_agree` column | Flags disagreement between GTEx and CellxGene for manual review |
| 5 | CPM back-transform + percentile column | Interpretable expression level; top 15% sorted first |
| 6 | F4 removed | Desirable biomarker profiles no longer dropped |
| 7 | TOP_N_FOR_EXTERNAL raised to 200 | RSPO3 and other validated candidates now reach external validation |
| 8 | Lesion-specificity distances as TOPSIS input | EuE no longer dominates via entropy weighting; direction-agnostic |
| 9 | UP/DOWN classification uses mean(EcP, EcO) | Fixes V1 bug where only EcO determined direction |
| 10 | EcPA excluded from TOPSIS criteria | Adjacent peritoneum is neither lesion nor healthy control |
| 11 | No penalty for genes outside TOP_N | Removes arbitrary 0.1x score multiplier from V1 |

---

## License

MIT License. See LICENSE file.

## Citation

If you use this pipeline, please cite the primary data source:

Tan, Yuliana et al. "Single-cell analysis of endometriosis reveals a coordinated transcriptional programme driving immunotolerance and angiogenesis across eutopic and ectopic tissues." *Nature Cell Biology* 24.8 (2022): 1306–1318. doi:10.1038/s41556-022-00961-5
