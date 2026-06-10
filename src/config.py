import os

# ==========================================
# V2 CONFIGURATION
# ==========================================

# 1. External Data Toggles
ENABLE_EXTERNAL_SPECIFICITY = True    # GTEx integration (was False in V1)
ENABLE_CELLXGENE = True               # CellxGene off-target integration
TOP_N_FOR_EXTERNAL = 200              # Raised from 50 in V1

# 2. Expression Thresholds
CPM_PERCENTILE_THRESHOLD = 0.85       # Top 15% by CPM = high-confidence expression

# 3. Dataset Configuration
CELL_TYPES_INFO = [
    (0,  'dS2'),
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
PLOTS_DIR  = os.path.join(OUTPUT_DIR, "plots")
RUNS_DIR   = os.path.join(OUTPUT_DIR, "runs")
CACHE_DIR  = ".api_cache"

# 5. GTEx Configuration
# Tissues excluded from off-target burden (reproductive/target-adjacent)
TARGET_TISSUES = ['Ovary', 'Uterus', 'Vagina', 'Fallopian Tube', 'Cervix Uteri']
GTEX_UTERUS_TISSUE_ID = 'Uterus'     # Reference tissue for log2-ratio normalization

# 6. CellxGene Configuration
# Data is fetched via the cellxgene-census Python API (requires Python < 3.13).
# Install: conda create -n endo_pipeline python=3.12 && pip install cellxgene-census
CELLXGENE_CENSUS_VERSION   = "2025-11-08"   # Pin to stable release (avoids version warning)
CELLXGENE_TIMEOUT_SECONDS  = 300            # Abort CellxGene fetch after 5 min if unresponsive
CELLXGENE_REPRODUCTIVE_TISSUES = [
    'uterus', 'ovary', 'vagina', 'fallopian tube', 'cervix', 'endometrium'
]
