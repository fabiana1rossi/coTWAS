# Configuration file for coTWAS analysis
# Paths and parameters

# Input paths
METADATA_PATH <- "./coTWAS/metadata.RData"
COMMON_GENES_PATH <- "./coTWAS/Common.Genes.Models_updated120924.RData"
UNIQUE_GENES_PATH <- "./coTWAS/Unique.Genes.Models_updated120924.RData"
GENE_INFO_PATH <- "./gene_annot.rds"

# Prediction paths
PREDICTIONS_BASE_PATH <- "./PGC3/predictions"
EPIXCAN_PATH <- file.path(PREDICTIONS_BASE_PATH, "EpiXcan")
CIS_PATH <- file.path(PREDICTIONS_BASE_PATH, "CIS")
MODULE_PATH <- file.path(PREDICTIONS_BASE_PATH, "MODULE/module_averaged")
INGENE_PATH <- file.path(PREDICTIONS_BASE_PATH, "INGENE/ingene_averaged")

# Output paths
OUTPUT_BASE_PATH <- "./coTWAS"
INPUT_PATH <- file.path(OUTPUT_BASE_PATH, "input")
OUTPUT_PATH <- file.path(OUTPUT_BASE_PATH, "output")

# Analysis parameters
TISSUES <- c("dlpfc", "hippo", "caudate", "amygdala", "sACC", "dACC")
MHC_REGION <- data.frame(
  chr = "chr6",
  start = 25726063,
  end = 33400644
)

# Create output directories if they don't exist
dir.create(INPUT_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_PATH, recursive = TRUE, showWarnings = FALSE) 