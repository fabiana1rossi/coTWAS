#####################################################################
# INGENE Training Pipeline
#####################################################################
# This script implements the INGENE (Inter-Network Gene Expression) approach
# for training gene expression prediction models.
#
# Key Features:
# - Uses elastic net regression to predict gene expression
# - Leverages gene co-expression networks for feature selection
# - Uses imputed expression of co-expression partners as predictors
# - Supports multiple brain regions and network modules
# - Performs cross-validation for model evaluation
# - Generates consolidated output files per network
# - Outputs comprehensive model summaries matching CIS/EpiXcan format
#
# Input Requirements:
# - Imputed expression data (EpiXcan and CIS predictions)
# - Network module definitions (co-expression networks)
# - Real expression data for training
# - Gene annotations
# - Covariates data (optional, for expression adjustment)
#
# Output Files:
# - model_summary.txt: Performance metrics for all genes in network (CIS/EpiXcan format)
# - model_weights.txt: Feature weights for prediction models (CIS/EpiXcan format)
# - prediction_results.txt: Predicted vs actual expression values
#
# Parameters:
# - alpha: Elastic net mixing parameter (0.5 = balanced between ridge and lasso)
# - n_folds: Number of cross-validation folds (default: 4)
# - n_train_test_folds: Number of train-test splits (default: 4)
# - include_cis: Whether to include cis-regulatory variants (not used in INGENE)
#####################################################################

# Required packages
required_packages = c(
    "purrr", "parallel", "tidyverse", "IRanges", 
    "limma", "glmnet"
)

# Install and load required packages
for(pkg in required_packages) {
    if(!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
        library(pkg, character.only = TRUE)
    }
}

## Set working directory
# Set working directory to project root
# setwd("path/to/project/root")

# Source helper functions
source("./STEP7.Training_INGENE/training_helper_functions.R")

# String concatenation operator
"%&%" = function(a,b) paste(a,b, sep='')

#################################################
## Configuration
#################################################

CONFIG = list(
    # Input directories and files
    input = list(
        # Directory containing imputed expression data
        imputed_data_dir  = "./predictions/",
        imputed_data_file = "CIS_Epi_Imputed.RData",
        
        # Network modules definition file
        modules_file = "./data/coexpression_network.rds", ## your coexpression networks
        
        # Training data location
        training_data = "./rnaseq/your_rnaseq.RData",
        
        # Cross validation folds file pattern
        cv_folds_pattern = "./data/cv_4fold_INGENE_%s_ids.RData" ## your fold ids
    ),
    
    # Output directory structure - Consolidated per network
    output = list(
        # Base output directory
        base_dir = "./STEP7.Training_INGENE/output/",
        
        # File patterns for consolidated network results
        files = list(
            # Summary file containing all model metrics
            network_summary = function(network_name) {
                paste0(network_name, "_summary.txt")
            },
            # Weights file containing all model coefficients
            network_weights = function(network_name) {
                paste0(network_name, "_weights.txt")
            },
            # Prediction results file
            network_predictions = function(network_name) {
                paste0(network_name, "_predictions.txt")
            }
        )
    ),
    
    # Brain regions to process
    brain_regions = c("dlpfc"),
    
    # Model training parameters
    model = list(
        alpha = 0.5,                # Elastic net mixing (0.5 = balanced L1/L2)
        n_folds = 4,               # Number of CV folds
        n_train_test_folds = 4,    # Number of train-test splits
        seed = 070891              # Random seed for reproducibility
    ),
    
    # Computation settings
    compute = list(
        n_cores = 16              # Number of CPU cores for parallel processing
    ),
    
    # Model filtering settings
    filters = list(
        include_cis = FALSE,       # Whether to include cis-regulatory variants
        custom_filter = NULL       # Additional custom filtering criteria
    )
)

# Function to validate configuration
validate_config = function(config) {
    # Check if required directories exist
    required_dirs = c(
        dirname(config$input$imputed_data_dir),
        dirname(config$input$modules_file),
        dirname(config$output$base_dir)
    )
    
    missing_dirs = required_dirs[!dir.exists(dirname(required_dirs))]
    if(length(missing_dirs) > 0) {
        stop("Missing required directories:\n", 
             paste(missing_dirs, collapse="\n"))
    }
    
    # Check if required files exist
    required_files = c(
        file.path(config$input$imputed_data_dir, config$input$imputed_data_file),
        config$input$modules_file
    )
    
    missing_files = required_files[!file.exists(required_files)]
    if(length(missing_files) > 0) {
        stop("Missing required files:\n",
             paste(missing_files, collapse="\n"))
    }
    
    # Validate computational parameters
    if(config$compute$n_cores > parallel::detectCores()) {
        warning("Requested cores (", config$compute$n_cores, 
                ") exceeds available cores (", parallel::detectCores(), ")")
    }
    
    return(TRUE)
}

# Validate configuration before running
validate_config(CONFIG)

#################################################
## Main Analysis Functions
#################################################

#' Process a single brain region's data for model training
#' @description Handles the complete training pipeline for one brain region:
#' 1. Creates output directories
#' 2. Loads and prepares training data
#' 3. Processes each network module
#' 4. Trains models for genes in each module
#' 5. Saves consolidated results per network
#' 
#' @param brain_region Brain region name 
#' @param data_imputed List containing imputed expression data for all regions
#' @param moduleList List of network module definitions
#' @param covs Covariates data frame
#' @param gene_annot Gene annotations
#' @param config Configuration parameters
#' @return NULL (saves results to files)
process_brain_region = function(brain_region, data_imputed, moduleList, covs, gene_annot, config) {
    message(sprintf("Processing brain region: %s", brain_region))
    
    # Create output directories
    dir.create(file.path(config$output$base_dir,"summary", brain_region), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(config$output$base_dir,"weights", brain_region), recursive = TRUE, showWarnings = FALSE)
    
    # Get training data
    train_real    = data_imputed[[brain_region]]$training.real
    train_imputed = data_imputed[[brain_region]]$training.predicted
 
    # Get network modules
    net_names = names(moduleList)
    
    # Process each network module
    results = map(net_names, function(net_name) {
        message(sprintf("Processing network: %s", net_name))
        
         
        # Define output files for this network
        network_files = list(
            summary = file.path(config$output$base_dir,"summary", brain_region,
                              config$output$files$network_summary(net_name)),
            weights = file.path(config$output$base_dir,"weights", brain_region,
                              config$output$files$network_weights(net_name))
        )
        if(file.exists(network_files[[1]]) && file.exists(network_files[[2]])){return(message(sprintf(" Network %s in %s already done :)!", net_name, brain_region)))}

        # Train models for this network using updated helper function
        INGENE.EnTraining(
            net_name      = net_name,
            brain.region  = brain_region,
            net.data      = moduleList,
            train.expr    = train_real,
            imputed.train = train_imputed,
            geneAnnot     = gene_annot,
            covs_df       = NULL,  # Pass covariates for adjustment
            output_files  = network_files,
            filter        = config$filters$custom_filter,
            include.cis   = config$filters$include_cis,
            no_cores      = config$compute$n_cores
        )
    })
}

#################################################
## Pipeline 
#################################################

#' Main pipeline execution function
#' @description Orchestrates the complete INGENE training pipeline:
#' 1. Initializes random seed and logging
#' 2. Loads and preprocesses input data
#' 3. Processes each brain region
#' 4. Records execution time and summary statistics
main = function() {
    message("Initializing INGENE training pipeline...")
    start_time = Sys.time()
    set.seed(CONFIG$model$seed)
    
    # Load data
    message("Loading input data...")
    load(file.path(CONFIG$input$imputed_data_dir, CONFIG$input$imputed_data_file))
    moduleList = readRDS(CONFIG$input$modules_file)
    
    # Filter module list (remove the grey module)
    moduleList = filter_module_list(moduleList)
    
    ## Get covs and gene annot
    training.data = get(load(CONFIG$input$training_data))
    
    # Process each brain region
    for(brain_region in CONFIG$brain_regions) {
        message(sprintf("\nProcessing brain region: %s", brain_region))
        
        ## Get covs and annot for region
        covs       = as.data.frame(training.data[[brain_region]]$colData)
        gene_annot = as.data.frame(training.data[[brain_region]]$rowData)
        
        ## Process region with updated function signature
        process_brain_region(brain_region, data_imputed = data.imputed, moduleList, covs, gene_annot, config = CONFIG)
    }
    
    # Report completion
    end_time = Sys.time()
    duration = difftime(end_time, start_time, units="mins")
    message(sprintf("\nPipeline completed in %.1f minutes", as.numeric(duration)))
}

# Run pipeline
main()


