#####################################################################
# INGENE Training Data Preparation Pipeline
#####################################################################
# This script prepares the input data for INGENE training by:
# 1. Loading and processing imputed expression data from different sources
# 2. Combining CIS and EpiXcan predictions
# 3. Creating training/testing datasets
#####################################################################

#################################################
## Setup and Configuration
#################################################

# Required packages
required_packages = c(
    "tidyverse", "limma", "IRanges", "purrr","data.table"
)

# Install and load required packages
for(pkg in required_packages) {
    if(!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
        library(pkg, character.only = TRUE)
    }
}

## Set working directory
# setwd("path/to/project/root")

# Source helper functions
source("./STEP7.Training_INGENE/prepare_training_data_helper_functions.R")

#################################################
## Configuration
#################################################

CONFIG = list(
    # Base directories
    base_dir = getwd(),
    
    # Input data paths
    input   = list(
        # Testing data
        testing_data = list(
            data = "./rnaseq/your_rnaseq.RData",
            epixcan = list(
                dir = "./predictions/EpiXcan",
                pattern = "EpiXcan_%s_testing_genotype_name_predicted.txt.gz",
                output = "EpiXcan_testing_genotype_name_IMPUTED.RData"
            ),
            cis = list(
                dir = "./predictions/CIS",
                pattern = "CIS_%s_testing_genotype_name_predicted.txt.gz",
                output = "CIS_testing_genotype_name_IMPUTED.RData"
            )
        ),
        
        # training data
        training_data = list(
            data = "./rnaseq/your_rnaseq.RData",
            epixcan = list(
              dir = "./predictions/EpiXcan",
              pattern = "EpiXcan_%s_testing_genotype_name_predicted.txt.gz",
              output = "EpiXcan_testing_genotype_name_IMPUTED.RData"
            ),
            cis = list(
              dir = "./predictions/CIS",
              pattern = "CIS_%s_testing_genotype_name_predicted.txt.gz",
              output = "CIS_testing_genotype_name_IMPUTED.RData"
            )
        )
    ),
    
    # Output paths
    output = list(
        dir = "./predictions/",
        file = "CIS_Epi_Imputed.RData"
    ),
    
    # Brain regions to process
    brain_regions_cis = c(
        "dlpfc","amygdala","hippo","acc"
    ),
    
    brain_regions_epixcan = c(
      "dlpfc","caudate","hippo","acc"
    )
)

#################################################
## Main Functions
#################################################

#' Load and process testing  data
#' @param config Configuration list
#' @return List of  expression data by tissue
load_testing_expr = function(config) {
    testing_data = get(load(config$input$testing_data$data))
    
    return(map(names(testing_data) %>% set_names(.,.), function(tissue) {
        return(testing_data[[tissue]]$assays$expression)}))
  
}

#' Load and process imputed expression data
#' @param config Configuration list
#' @param data_type Type of data to load ("epixcan" or "cis")
#' @param source Data source ("testing" or "traininig")
#' @return List of imputed expression data by brain region
load_imputed_data = function(config, data_type, source, model) {
    input_config  = config$input[[source]][[data_type]]
    
    if(model=="EpiXcan"){brain_regions = config$brain_regions_epixcan}
    if(model=="CIS"){brain_regions = config$brain_regions_cis}
    
    predicted_lst = lapply(brain_regions, function(region) {
        file_path = file.path(input_config$dir, sprintf(input_config$pattern, region))
        tmp       = fread(file_path)
        return(as.data.frame(tmp))
   })
    names(predicted_lst) = brain_regions
    
    output_path = file.path(dirname(input_config$dir), input_config$output)
    save(predicted_lst, file = output_path)
    
    return(predicted_lst)
}

#' Process data for each brain region
#' @param config Configuration list
#' @param testing_expr testing expression data
#' @param epixcan_testing EpiXcan testing predictions
#' @param cis_testing CIS testing predictions
#' @param epixcan_traininig EpiXcan traininig predictions
#' @param cis_traininig CIS traininig predictions
#' @param traininig_data traininig expression data
#' @return List of processed data by brain region
process_brain_regions = function(config, testing_data, epixcan_testing, cis_testing,
                               epixcan_traininig, cis_traininig, training_data) {
    data.imputed = map(config$brain_regions_cis %>% set_names(.,.), function(region) {
        message(sprintf("Processing region: %s", region))
        
        # Get data for current region
        data = prepare_region_data(
            testing.real.data              = testing_data[[region]],
            epixcan_testing.imputed.data   = epixcan_testing[[region]],
            cis_testing.imputed.data       = cis_testing[[region]],
            epixcan_training.imputed.data  = epixcan_traininig[[region]],
            cis_training.imputed.data      = cis_traininig[[region]],
            training.real                  = training_data[[region]]
        )
        
        return(data)
    })
    
    return(data.imputed)
}

#################################################
## Run Pipeline
#################################################

main = function() {
    message("Starting data preparation pipeline...")
    start_time = Sys.time()
    
    # Load testing expression data &  training expression data
    training_data = get(load(CONFIG$input$training_data$data))
    
    ## Extract expression
    training_data = map(names(training_data) %>% set_names(.,.), function(region){
      tmp = as.data.frame(training_data[[region]][["assays"]][["expression"]])
      
      return(tmp)
    })
    
    testing_data=training_data
    
    # Load imputed data
    epixcan_testing      = epixcan_traininig = load_imputed_data(CONFIG, "epixcan", "testing_data", "EpiXcan")
    cis_testing          = cis_traininig     = load_imputed_data(CONFIG, "cis", "testing_data", "CIS")
    
    # epixcan_traininig    = load_imputed_data(CONFIG, "epixcan", "training_data")
    # cis_traininig        = load_imputed_data(CONFIG, "cis", "training_data")
    
    # Process each brain region
    data.imputed = process_brain_regions(
        CONFIG, testing_data, epixcan_testing, cis_testing,
        epixcan_traininig, cis_traininig, training_data
    )
    
    # Save final output
    output_path = file.path(CONFIG$output$dir, CONFIG$output$file)
    save(data.imputed, file = output_path)
    
    end_time = Sys.time()
    message(sprintf("Pipeline completed in %s", end_time - start_time))
}

# Run pipeline
main() 
