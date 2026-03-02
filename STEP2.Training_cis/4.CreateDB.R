#####################################################################
# Create Prediction Database Pipeline
#####################################################################
# This script creates SQLite databases from the training results
# It performs the following steps:
# 1. Processes summary and weights files for each region
# 2. Creates SQLite databases for prediction models
# 3. Saves databases in specified output directory
#####################################################################

#################################################
## Setup and Configuration
#################################################

# Load required libraries
required_packages = c(
  "dplyr", "purrr", "RSQLite", "magrittr", "limma",
  "data.table", "DBI", "parallel"
)

# Install and load packages
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


## set working directory
# setwd("")

# Source helper functions
source("./STEP2.Training_cis/database_helper_functions.R")

#################################################
## Configuration
#################################################
base_dir = getwd()

CONFIG = list(
  summary_dir       = file.path(base_dir, "STEP2.Training_cis/output_EpiXcan/summary/"),
  weights_dir       = file.path(base_dir, "STEP2.Training_cis/output_EpiXcan/weights/"),
  db_dir            = file.path(base_dir, "STEP2.Training_cis/output_EpiXcan/database/"),
  brain_regions     = c('dlpfc'),
  dataset_training  = 'your_dataset_training_name',
  model             = 'EpiXcan',
  n_cores           = 1
)

# Create output directories if they don't exist
if (!dir.exists(CONFIG$db_dir)) dir.create(CONFIG$db_dir, recursive = TRUE)

#################################################
## Run Database Creation Pipeline
#################################################

# Process each brain region
for(brain_region in CONFIG$brain_regions) {
  message("Processing brain region: ", brain_region)
  
 
  # Create databases for each network
  db_files = makeModelDB1(
    brain.region       = brain_region,
    file_input_summ    = CONFIG$summary_dir,
    file_input_weights = CONFIG$weights_dir,
    drv_file_path      = CONFIG$db_dir,
    dataset_training   = CONFIG$dataset_training,
    model              = CONFIG$model
  )
  
  message("Created database files for brain region ", brain_region, ":")
  print(db_files)
}

message("Database creation pipeline completed successfully!") 
