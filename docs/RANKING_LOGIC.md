# Biomarker Pipeline Ranking: TOPSIS Algorithm

This document explains the mathematical logic behind the ranking strategy used in the pipeline, which you can use for your presentation.

## 1. What is TOPSIS?
**TOPSIS** stands for *Technique for Order of Preference by Similarity to Ideal Solution*. It is a multi-criteria decision-analysis method that ranks candidates based on their geometric distance to a theoretical "Ideal" and "Anti-Ideal" solution.

## 2. Why use TOPSIS for Biomarkers?
In this pipeline, we have conflicting goals:
1.  **Maximize Signal:** High absolute expression change in Ectopic lesions (EcP, EcO).
2.  **Minimize Signal:** Expression change as close to **0** as possible in Eutopic tissue (EuE) to ensure specificity.

TOPSIS allows us to balance these goals mathematically without relying on simple averages, which can be easily skewed by single extreme values.

## 3. The Mathematical Workflow

### Step 1: Defining the Targets
We define two imaginary genes in our multi-dimensional space:
*   **The Ideal Biomarker ($V^+$):** Has the maximum observed logFC in lesions across the entire dataset and **exactly 0** logFC in eutopic tissue.
*   **The Anti-Ideal Biomarker ($V^-$):** Has **0** signal in lesions and the maximum observed logFC in eutopic tissue (the worst possible specificity).

### Step 2: Vector Normalization
Since different tissue comparisons might have different ranges of values, we normalize the data so that every comparison is on the same scale (0 to 1). This prevents one tissue group from dominating the score just because it has higher variance.

$$n_{ij} = \frac{x_{ij}}{\sqrt{\sum x_{ij}^2}}$$

### Step 3: Automated Weighting (Entropy Weighting Method)
Instead of manual weights, the pipeline now automatically calculates the "Ideal Weights" based on the data's inherent diversification.
1.  **Information Entropy:** We calculate the entropy for each tissue comparison. High entropy means the values are very similar (low information).
2.  **Diversification Degree:** $d_j = 1 - e_j$. This measures how well a specific tissue comparison can "discriminate" between different genes.
3.  **Objective Weights:** The final weight $w_j$ is the normalized diversification degree. A column that shows high variance and clear separation between candidates is automatically given a higher weight.

This ensures that the ranking is purely data-driven and prioritizes the most informative biomarkers.

### Step 4: Geometric Distance Calculation
For every gene, we calculate:
*   **$S^+$:** The Euclidean distance to the **Ideal Solution**.
*   **$S^-$:** The Euclidean distance to the **Anti-Ideal Solution**.

### Step 5: Final Score (Relative Closeness)
The final **Score** is calculated as:

$$C_i = \frac{S^-_i}{S^+_i + S^-_i}$$

*   A score of **1.0** means the gene is exactly at the Ideal point.
*   A score of **0.0** means the gene is exactly at the Anti-Ideal point.

## 4. Interpretation
Genes at the top of the `Rank` (highest Score) are those that most effectively "stretch" the distance between being perfect and being terrible. This provides a robust, reproducible, and mathematically defensible way to select the best biomarker candidates for further validation.
