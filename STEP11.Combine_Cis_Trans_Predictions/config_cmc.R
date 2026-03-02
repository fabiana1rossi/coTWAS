# Script: config.R
# Purpose: Centralizes all configuration parameters for the linear modeling analysis pipeline

###########################
# Base Directory Settings #
###########################

# Set working directory and base path
# BASE_DIR should be set to the project root directory
# BASE_DIR = "path/to/project/root"
BASE_DIR = getwd()

###########################
# Data Paths             #
###########################

# Input data paths
DATA_PATHS = list(
  # RNA-seq data
  rnaseq = list(
    testing = file.path(BASE_DIR, "rnaseq/cmc_rnaseq.RData")
  ),
  
  # Predictions directory
  predictions = list(
    base = file.path(BASE_DIR, "predictions/LIBD_training/CMC_EUR/"),
    combined = file.path(BASE_DIR, "predictions/LIBD_training/CMC_EUR/Combined")
  )
)

###########################
# Analysis Parameters    #
###########################

# Brain regions to analyze
BRAIN_REGIONS = c("dlpfc","dACC","sACC")

# Dataset name
DATASET_NAME = "CMC"

# Threshold setting
THRESHOLD = "yes"

###########################
# File Naming Patterns   #
###########################

FILE_PATTERNS = list(
  # MODULE files
  module = list(
    predictions = "{brain_region}_AveragedNetworks_predicted.txt",
    performance = "MODULE_replicable_{dataset}_performance_{brain_region}.RData"
  ),
  
  # EpiXcan files
  epixcan = list(
    predictions = "EpiXcan_{brain_region}_{dataset}_predicted.txt.gz",
    performance = "EpiXcan_{brain_region}_{dataset}_performance_selected.RData"
  ),
  
  # INGENE files
  ingene = list(
    predictions = "{brain_region}_AveragedNetworks_predicted.txt",
    performance = "INGENE_replicable_{dataset}_performance_{brain_region}.RData" 
  ),
  
  # CIS files
  cis = list(
    predictions = "CIS_{brain_region}_{dataset}_predicted.txt.gz",
    performance = "CIS_{brain_region}_{dataset}_performance_selected.RData"
  ),
  
  # Output files
  output = list(
    predictions_lst = "predictions_testing_data_lst.RData",
    model_results = "Common_Predictions_GTEx_{brain_region}_replicable.RData"
  )
)

###########################
# Helper Functions       #
###########################

#' Get file path with proper formatting
#' @param pattern File pattern from FILE_PATTERNS
#' @param ... Named parameters to substitute in pattern
#' @return Formatted file path
get_file_path = function(pattern, ...) {
  params = list(...)
  path = pattern
  for (param in names(params)) {
    path = gsub(paste0("\\{", param, "\\}"), params[[param]], path)
  }
  return(path)
}

#' Get full path for a specific file type
#' @param method Method name (module, epixcan, ingene, cis)
#' @param file_type Type of file (predictions, performance)
#' @param ... Additional parameters for path formatting
#' @return Full file path
get_method_file_path = function(method, file_type, ...) {
  pattern = FILE_PATTERNS[[method]][[file_type]]
  return(get_file_path(pattern, ...))
} 

