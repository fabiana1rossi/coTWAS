#####################################################################
# Brain Tissue Module coeQTL Analysis Pipeline
#####################################################################
# This script performs co-expression QTL analysis for brain tissue modules.
# It processes the output from the mapping pipeline to:
# 1. Perform cross-validation for coeQTL detection
# 2. Rank and prune coeQTLs based on significance
# 3. Save results for each network
#####################################################################



################################################
## Setup and Configuration
#################################################

# Load required libraries
required_packages = c(
    "purrr", "parallel", "limma", "pbapply", "robust", 
    "robustbase", "MASS", "sfsmisc", "magrittr", "dplyr", 
    "matrixStats"
)

# Install and load packages
for(pkg in required_packages) {
    if(!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
        library(pkg, character.only = TRUE)
    }
}

## Set working directory
# Set working directory to project root
# setwd("path/to/project/root")


#################################################
## Helper Functions  ##
#################################################
source("./STEP3.Training_MODULE/coEQTLs_helper_functions.R")


#####################
# Configuration
#################
base_dir = getwd()

CONFIG = list(
    rnaseq              = paste0(base_dir, "/rnaseq/your_rnaseq.RData"), ## Rdata expression & covariates file
    data_dir            =  paste0(base_dir, "/data/"),
    n_cores             = 1,  ## specify the number of cores to run the script in parallel (ubuntu machine)
    n_folds             = 4,  ## cross-validation folds
    brain_regions       = c('dlpfc'), ## Specify available regions
    flank_size          = 250000, ## windowd for LDpruning
    r_squared_threshold = 0.1     ## r2 for LDpruning
)

# Set paths
CONFIG$network_res_dir = CONFIG$output_dir = file.path(base_dir, "STEP3.Training_MODULE/output/") ## The module/network specific RData created in the previous step
                                                                                                  ## will be updated with coEQTLs results



#################################################
## Usage 
#################################################

run_coeqtl_pipeline(
        expr_data_covs = CONFIG$rnaseq,
        network_res_dir = CONFIG$network_res_dir,
        data_dir = CONFIG$data_dir,
        brain_regions = CONFIG$brain_regions,
        n_cores = CONFIG$n_cores,
        n_folds = CONFIG$n_folds,
        flank_size = CONFIG$flank_size,
        r_squared_threshold = CONFIG$r_squared_threshold
    )

