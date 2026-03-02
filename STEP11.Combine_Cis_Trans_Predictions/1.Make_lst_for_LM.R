# Script: 1.Make_lst_for_LM.R
# Purpose: Prepares testing data for linear modeling analysis by:
#   1. Loading and processing RNA-seq expression data
#   2. Processing predictions from different methods (MODULE, INGENE, EpiXcan, CIS)
#   3. Combining predictions and real data for each brain region


###########################
# Load Configuration     #
###########################
# Set working directory to load config
# setwd("path/to/project/root")
# Load configuration
source("./STEP11.Combine_Cis_Trans_Predictions/config.R")

# Set working directory (BASE_DIR is set in config.R)
setwd(BASE_DIR)

###########################
# Load Required Libraries #
###########################

# Data manipulation libraries
library(data.table)  
library(tidyverse)   
library(limma)       

# Load custom functions
source("./STEP11.Combine_cis_trans_predictions/LM_models_helper_functions.R")

###########################
# Load and Process Data  #
###########################

#' Load and Process Testing Data
#' @param data_path Path to the testing RNA seq data
#' @param tissues Vector of tissue names to analyze
#' @return List of expression data for each tissue
load_testing_data = function(data_path = DATA_PATHS$rnaseq$testing, 
                           tissues = BRAIN_REGIONS) {
  testing_data = get(load(data_path))
  testing_data$sACC = testing_data$dACC = testing_data$acc
  map(tissues %>% set_names(.,.), ~ {
    tissue = .x
    testing_data[[tissue]]$assays$expression
  })
}

# Load testing data expression data
testing_data = load_testing_data()

###########################
# Process Predictions    #
###########################

# Process testing data predictions for each brain region
testing_data_lst = map(names(testing_data), function(brain.region){
  print(paste("Processing brain region:", brain.region))
  
  # Prepare data for LM analysis
  lst = getDataforLM(
    dataset           = DATASET_NAME,  
    input.dir         = DATA_PATHS$predictions$base,
    true_expr         = testing_data[[brain.region]],
    threshold         = THRESHOLD,
    brain.region      = brain.region
  )
  ## Add real data 
  lst[["real.expr"]]  = testing_data[[brain.region]]
  return(lst)
})

# Name the results and save
names(testing_data_lst) = names(testing_data)
if(!dir.exists(file.path(DATA_PATHS$predictions$combined))){dir.create(file.path(DATA_PATHS$predictions$combined), recursive = T)}
save(testing_data_lst, 
     file = file.path(DATA_PATHS$predictions$combined, 
                     FILE_PATTERNS$output$predictions_lst))

print('Data processing completed successfully')


