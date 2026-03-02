#####################################################################
# INGENE testing data Prediction Pipeline - Network Averaging
#####################################################################
# This script implements the final step of the INGENE pipeline where
# predictions from multiple networks are combined. The process involves:
# 1. Loading predictions from different networks
# 2. Validating predictions against true expression data
# 3. Averaging predictions for genes present in multiple networks
# 4. Preserving unique predictions from single networks
# 5. Saving both averaged predictions and performance metrics
#####################################################################
## Set working directory
# Set working directory to project root
# setwd("path/to/project/root")

# Function to check if required files and directories exist
check_environment = function(config) {
  # Check if predictions directory exists
  if (!dir.exists(config$input$predictions_dir)) {
    stop("Predictions directory not found: ", config$input$predictions_dir)
  }
  
  # Check if each required subdirectory exists
  for (cis_type in config$input$cis_types) {
    for (region in config$input$regions) {
      dir_path = file.path(config$input$predictions_dir, cis_type, region)
      if (!dir.exists(dir_path)) {
        warning("Directory not found and will be skipped: ", dir_path)
      }
    }
  }
}

# Function to safely read prediction files
safe_read_predictions = function(file_path) {
  tryCatch({
    data = fread(file_path)
    data = as.data.frame(data)
    if (nrow(data) == 0) {
      warning("File is empty: ", file_path)
      return(NULL)
    }
    if (ncol(data) < 2) {  # At minimum should have IID, FID, and one gene
      warning("File has insufficient columns: ", file_path)
      return(NULL)
    }
    return(data)
  }, error = function(e) {
    warning("Error reading file ", file_path, ": ", e$message)
    return(NULL)
  })
}

base_dir = getwd()

# Load helper functions for averaging and validation
required_sources = c(
  './STEP10. Average_Network_Predictions/average_net_helper_functions.R',
  './STEP8.Validate_cis_predictions_testing_data/validate_predictions_helper_functions.R'
)

for (source_file in required_sources) {
  if (!file.exists(source_file)) {
    stop("Required source file not found: ", source_file)
  }
  source(source_file)
}

# Required packages for data manipulation and analysis
required_packages = c(
  "dplyr", "purrr","limma","data.table"
)

# Install and load required packages
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    message("Installing package: ", pkg)
    install.packages(pkg)
    if(!require(pkg, character.only = TRUE)) {
      stop("Failed to install and load package: ", pkg)
    }
  }
}

#################################################
## Configuration Settings
#################################################

CONFIG = list(
  input = list(
    # Data types to process
    cis_types = c("EpiXcan", "CIS"),  # Types of predictions to average
    regions = c('dlpfc'),             # Brain regions to analyze
    
    # Directory containing INGENE predictions
    predictions_dir = './predictions/INGENE/',
    
    ## cores
    n_cores = 8
  ),
  
  output = list(
    # Output directory for averaged predictions
    base_dir = './predictions/INGENE/INGENE_averaged/',
    file_prefix = "INGENE_averaged_"
  )
)

# Check environment before proceeding
check_environment(CONFIG)

# Create output directory structure
for(region in CONFIG$input$regions) {
  for(cis in CONFIG$input$cis_types) {
    dir_path = file.path(CONFIG$output$base_dir, cis, region)
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dir_path)) {
      stop("Failed to create output directory: ", dir_path)
    }
  }
}

#' Load and Process Testing Data
#' @param data_path Path to the testing RNA seq data
#' @param tissues Vector of tissue names to analyze
#' @return List of expression data for each tissue
load_testing_data = function(data_path = './rnaseq/your_testing_rnaseq.RData', 
                             tissues = c("dlpfc")) {
  testing_data = get(load(data_path))
  map(tissues %>% set_names(.,.), ~ {
    tissue = .x
    testing_data[[tissue]][['assays']][['expression']]
  })
}

#################################################
## Main Pipeline Function
#################################################

main = function() {
  start_time = Sys.time()
  success_count = 0
  error_count = 0
  
  # Process each prediction type (CIS/EpiXcan)
  for(cis_type in CONFIG$input$cis_types) {
    message("\nProcessing ", cis_type)
    
    
    # Process each brain region
    for(region in CONFIG$input$regions) {
      message("\nProcessing region: ", region)
      
      if(file.exists(paste0(file.path(CONFIG$output$base_dir, cis_type), "/", region, "_",cis_type,"_AveragedNetworks_predicted.txt"))){message("ALREDEY DONE");next}
      
      tryCatch({
        # Get prediction files for current region
        pred_dir = file.path(CONFIG$input$predictions_dir, cis_type, region)
        if (!dir.exists(pred_dir)) {
          warning("Directory not found: ", pred_dir)
          error_count = error_count + 1
          next
        }
        
        pred_files = list.files(pred_dir, pattern = "_predictions.txt$", full.names = TRUE)
        
        if(length(pred_files) == 0) {
          warning("No prediction files found for ", region, " in ", cis_type)
          error_count = error_count + 1
          next
        }
        
        message("Found ", length(pred_files), " prediction files")
        
        # Read all prediction files
        predictions_list = lapply(pred_files, safe_read_predictions)
        predictions_list = predictions_list[!sapply(predictions_list, is.null)]
        
        if(length(predictions_list) == 0) {
          warning("No valid prediction data found for ", region, " in ", cis_type)
          error_count = error_count + 1
          next
        }
        
        message("Successfully loaded ", length(predictions_list), " prediction files")
        
        # Load testing data
        testing_data_true_expr = load_testing_data()
        #testing_data_true_expr$sACC = testing_data_true_expr$dACC = testing_data_true_expr$acc
        
        ## get network names
        net.names = gsub("_predictions.txt","", basename(pred_files))
        net.names = gsub(paste0("INGENE_",region,"_"),"", net.names)
        
        
        # Average predictions across networks
        averaged_predictions = ensembleAcrossNet(
          model        = cis_type,
          brain.region = region,
          net.names    = net.names,
          true.expr    = testing_data_true_expr[[region]], 
          outputdir    = file.path(CONFIG$output$base_dir, cis_type),
          outname      = "predicted",
          network_predictions = predictions_list, 
          testing_geno = "LIBD",
          pipeline     = "INGENE",
          genestokeep  = NULL,
          n_cores      = CONFIG$input$n_cores
        )
        
        rm(predictions_list);gc()
        
        if (is.null(averaged_predictions)) {
          warning("Failed to generate averaged predictions for ", region, " in ", cis_type)
          error_count = error_count + 1
          next
        }
        
        success_count = success_count + 1
        message("Successfully processed ", region, " in ", cis_type)
        
      }, error = function(e) {
        error_count = error_count + 1
        warning("Error processing ", region, " in ", cis_type, ": ", e$message)
      })
    }
  }
  
  end_time = Sys.time()
  message("\nPipeline Summary:")
  message("----------------")
  message("Total time: ", difftime(end_time, start_time, units="mins"), " minutes")
  message("Successful processes: ", success_count)
  message("Failed processes: ", error_count)
  if (error_count > 0) {
    warning("\nSome processes failed. Check the warnings above for details.")
  }
}

# Execute the pipeline
tryCatch({
  main()
}, error = function(e) {
  message("\nCritical error in pipeline execution:")
  message(e$message)
  message("\nStack trace:")
  print(sys.calls())
}) 
