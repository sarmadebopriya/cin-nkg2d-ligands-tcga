"""Central settings. Edit this file, not the scripts.

Two environment variables let you point the pipeline elsewhere without editing code:
  CIN_DATA_DIR  folder holding raw/ and cin70_genes.txt   (default: ./data)
  CIN_OUT_DIR   folder for processed tables, results, figures (default: ./)
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_DIR = Path(os.environ.get("CIN_DATA_DIR", ROOT / "data"))
OUT_DIR = Path(os.environ.get("CIN_OUT_DIR", ROOT))

RAW_DIR = DATA_DIR / "raw"
PROC_DIR = OUT_DIR / "data_processed"
RESULTS_DIR = OUT_DIR / "results"
FIG_DIR = OUT_DIR / "figures"

# ---------------------------------------------------------------- input files
# Expression: tab-separated, first column = gene symbol, remaining columns = samples
# (UCSC Xena orientation). Values may be log2(x+1) already; the loader checks.
# Default is the PanCanAtlas batch-corrected file from the GDC. A Xena download saved as
# expression.tsv.gz also works: the loader uses whichever of these exists.
EXPR_FILE_CANDIDATES = [
    RAW_DIR / "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",
    RAW_DIR / "expression.tsv.gz",
    RAW_DIR / "expression.tsv",
]
EXPR_FILE = next((p for p in EXPR_FILE_CANDIDATES if p.exists()), EXPR_FILE_CANDIDATES[0])

# Source files for 00_build_phenotype.py (all from the GDC PanCanAtlas pages)
ANEUPLOIDY_FILE = RAW_DIR / "PANCAN_ArmCallsAndAneuploidyScore_092817.txt"
PURITY_FILE = RAW_DIR / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
MAF_FILE = RAW_DIR / "mc3.v0.2.8.PUBLIC.maf.gz"          # optional, for TP53 status
PAM50_FILE = RAW_DIR / "brca_clinical.tsv"               # optional, for the BRCA basal-like subset

PANCAN_CLINICAL_FILE = RAW_DIR / "pancan_clinical.tsv"   # optional: cBioPortal pan-cancer clinical table (FGA, MSI)
MSI_HIGH_MANTIS = 0.4     # MANTIS score above this is treated as MSI-high (commonly used cutoff); robustness only

# Fallback if downloads stall: set to ["BRCA"] to run the whole pipeline on breast cancer only.
ONLY_TYPES = None

# Phenotype: tab-separated, one row per sample (or per patient, see PHENO_KEY).
PHENO_FILE = RAW_DIR / "phenotype.tsv"
PHENO_KEY = "sample"          # "sample" if rows are TCGA sample barcodes, "patient" if 12-character patient IDs

# OPTIONAL. Carter 2006 CIN70 list, one gene symbol per line. If the file is absent, the CIN70 check is skipped.
CIN70_FILE = DATA_DIR / "cin70_genes.txt"

# Column names in the phenotype file. 00_build_phenotype.py writes exactly these names, so you
# only need to edit this if you supply your own phenotype file.
# Leave a value as None if you do not have that column (tp53_mut and pam50 are optional).
COLS = {
    "id": "sample",              # sample barcode or patient ID, matching PHENO_KEY
    "cancer_type": "cancer_type",
    "purity": "purity",
    "aneuploidy": "aneuploidy_score",
    "fga": "fga",                # fraction of genome altered
    "tp53_mut": "tp53_mut",      # 1 = mutated, 0 = wild type
    "pam50": "pam50",            # used only for the BRCA basal-like subset
    "msi": "msi_mantis",         # MSI MANTIS score, used only for the MSI-high exclusion check
}

# ---------------------------------------------------------------- gene lists
LIGANDS = ["MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "RAET1E", "RAET1G"]  # RAET1E = ULBP4, RAET1G = ULBP5
SENSING_EVASION = ["CGAS", "STING1", "TREX1", "ADAM10", "ADAM17"]
NK_MARKERS = ["NCR1", "KLRD1", "PRF1", "GZMB"]

# Older annotations use other symbols. The first name present in the expression file is used.
ALIASES = {
    "STING1": ["STING1", "TMEM173"],
    "RAET1E": ["RAET1E", "ULBP4"],
    "RAET1G": ["RAET1G", "ULBP5"],
    "CGAS": ["CGAS", "MB21D1", "C6ORF150"],
}
# Some older files label a gene "?" or an outdated symbol. Entrez IDs (the number after "|" in
# the row label) never change, so they are used first for these genes.
ENTREZ_TO_GENE = {"115004": "CGAS", "340061": "STING1"}
# The run stops if any of these cannot be found. Others only trigger a warning.
REQUIRED_GENES = ["MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "CGAS", "STING1"]

# ---------------------------------------------------------------- analysis choices
# Fixed before looking at results. Do not change them after seeing output.
PRIMARY_INSTABILITY = "aneuploidy"       # secondary: fga, tertiary: cin70
SECONDARY_INSTABILITY = ["fga", "cin70"]
PRIMARY_OUTCOME = "ligand_score"         # per-gene results are supplementary
MIN_N = 150                              # minimum tumors per cancer type after dropping missing values
HIGH_FRACTION = 1 / 3                    # "instability-high" = top third within each cancer type
MIN_GROUP = 20                           # minimum tumors per ligand-low / ligand-high group in Figure 2
MIN_N_SUBSET = 100                       # looser floor for the labeled BRCA basal-like subset
MIN_GROUP_SUBSET = 15
FDR_ALPHA = 0.05
PURITY_HIGH_QUANTILE = 0.5               # sensitivity analysis keeps tumors at or above this purity quantile
SEED = 20261124

# SKCM is mostly metastatic: prefer sample type 06 there, 01 everywhere else.
SAMPLE_TYPE_PRIORITY = {"default": ["01"], "SKCM": ["06", "01"]}

# Labeled subset analysed alongside the cancer types (not part of the FDR family).
BRCA_TYPE_LABEL = "BRCA"
BRCA_BASAL_VALUES = ["Basal", "Basal-like", "basal", "BRCA.Basal", "BRCA_Basal"]
