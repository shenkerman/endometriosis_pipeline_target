import os

# ==========================================
# CONFIGURATION
# ==========================================

# 1. External Data Toggles
ENABLE_EXTERNAL_SPECIFICITY = False  # Toggle GTEx integration
TOP_N_FOR_EXTERNAL = 50             # Limit external analysis to top candidates for speed

# 2. Ranking Weights (Used only as initial baseline or if automated weighting is off)
# Priority: logFC near 0 in EuE, high absolute logFC in lesions
RANKING_WEIGHTS = {
    'logFC.eueVSctrl': -2.0,
    'logFC.ecpVSctrl': 1.0,
    'logFC.ecoVSctrl': 1.0,
    'logFC.ecpaVSctrl': 0.5,
    'off_target_burden': -1.5  # Weight for GTEx-derived specificy (minimized)
}

# 3. Dataset Configuration
CELL_TYPES_INFO = [
    (0, 'dS2'), 
    (10, 'Prv-CCL19'), 
    (20, 'EC-tip'), 
    (30, 'EC-aPCV'), 
    (40, 'EC-PCV'), 
    (50, 'Mɸ1-LYVE1'), 
    (60, 'Mɸ4-infiltrated'), 
    (70, 'cDC2'), 
    (80, 'Treg'), 
    (90, 'B cell')
]

# 4. Directory Structure
OUTPUT_DIR = "results"
TABLES_DIR = os.path.join(OUTPUT_DIR, "tables")
PLOTS_DIR = os.path.join(OUTPUT_DIR, "plots")
RUNS_DIR = os.path.join(OUTPUT_DIR, "runs")
CACHE_DIR = ".api_cache"

# 5. GTEx Specific Configuration
# Tissues to EXCLUDE from "Off-target Burden" calculation (considered target-adjacent or ovary)
TARGET_TISSUES = ['Ovary', 'Uterus', 'Vagina', 'Fallopian Tube', 'Cervix Uteri']
