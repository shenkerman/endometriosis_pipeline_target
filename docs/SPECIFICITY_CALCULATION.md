# Off-target Burden Calculation: Technical Specification

This document provides a detailed breakdown of how the **Off-target Burden** metric is calculated and integrated into the biomarker pipeline.

## 1. Biological Objective
The goal of this metric is to quantify the "systemic noise" of a potential biomarker. An ideal biomarker should be nearly exclusive to the disease site (endometriotic lesions) or the tissue of origin (reproductive system). If a gene is highly expressed in other vital organs (e.g., liver, brain, heart), it is considered to have a high "off-target burden," making it a poor candidate for diagnostic tests or targeted therapies.

---

## 2. Data Sourcing (GTEx v8)
The pipeline retrieves real-world healthy tissue expression data from the **GTEx Portal (Genotype-Tissue Expression)**.

1.  **Identifier Mapping:** The pipeline maps HGNC Gene Symbols to **Gencode/Ensembl IDs** (e.g., `EFEMP1` $\rightarrow$ `ENSG00000115380.19`) using the GTEx Reference API. This ensures unambiguous data retrieval.
2.  **Median Expression:** We fetch the **Median TPM** (Transcripts Per Million) for the gene across all 54 tissue sites available in the GTEx v8 release.

---

## 3. The Calculation Formula
The Off-target Burden ($B$) is calculated as the sum of median expression across all healthy tissues, **excluding** the reproductive/target tissues defined in the configuration.

$$B_i = \sum_{t \in T_{off}} TPM_{it}$$

Where:
*   $i$ is the gene of interest.
*   $T_{off}$ is the set of all GTEx tissues **minus** the target tissues.
*   The default **Excluded Tissues** ($T_{target}$) are: `Ovary`, `Uterus`, `Vagina`, `Fallopian Tube`, and `Cervix Uteri`.

### Why Summation?
We sum the TPM values because a gene that is expressed at 10 TPM in 10 different organs is twice as "noisy" (systemically prevalent) as a gene expressed at 10 TPM in only 5 organs.

---

## 4. Integration into Ranking (TOPSIS)
Because TPM values (which can range from 0 to 100,000) are on a different scale than LogFC values (typically -5 to 5), the pipeline performs two normalization steps before ranking:

1.  **Vector Normalization:** The burden values are normalized to a 0-1 scale alongside the LogFC values.
2.  **Entropy Weighting:** The algorithm calculates the "Diversification Degree" of the burden column. If the candidates show very different levels of systemic noise, the algorithm automatically increases the weight of this metric in the final score.
3.  **Ideal Target:** In the TOPSIS model, the "Ideal Solution" for the Off-target Burden is defined as **0**.

---

## 5. Visual Interpretation
In the generated heatmaps (`all_cell_types_top_10_heatmap.png`), this value is shown in the **Magma** color spectrum:
*   **Deep Purple/Black (near 0):** High specificity (Ideal).
*   **Bright Orange/Yellow (High TPM):** Low specificity (Potential systemic noise).

This dual-validation approach ensures that the top-ranked candidates in the pipeline are not just statistically significant in the study, but biologically specific in the context of the whole human body.
