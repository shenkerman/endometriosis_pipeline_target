import os

# ==========================================
# CONFIGURATION
# ==========================================

# Priority: logFC near 0 in EuE, high absolute logFC in lesions
RANKING_WEIGHTS = {
    'logFC.eueVSctrl': -2.0,  # Negative weight to minimize absolute value
    'logFC.ecpVSctrl': 1.0,   # Positive weight to maximize absolute value
    'logFC.ecoVSctrl': 1.0,   # Positive weight to maximize absolute value
    'logFC.ecpaVSctrl': 1.0   # Positive weight to maximize absolute value
}

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

OUTPUT_DIR = "results"
TABLES_DIR = os.path.join(OUTPUT_DIR, "tables")
PLOTS_DIR = os.path.join(OUTPUT_DIR, "plots")
