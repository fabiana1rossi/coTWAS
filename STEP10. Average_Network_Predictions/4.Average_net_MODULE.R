#####################################################################
# MODULE testing data Prediction Pipeline - Network Averaging
#####################################################################
# This script implements the final step of the MODULE pipeline where
# predictions from multiple networks are combined. The process involves:
# 1. Loading predictions from different networks
# 2. Validating predictions against true expression data
# 3. Averaging predictions for genes present in multiple networks
# 4. Preserving unique predictions from single networks
# 5. Saving both averaged predictions and performance metrics
#####################################################################

# Function to check if required files and directories exist
check_environment = function(config) {
    # Check if predictions directory exists
    if (!dir.exists(config$input$predictions_dir)) {
        stop("Predictions directory not found: ", config$input$predictions_dir)
    }
    
}

# Function to detect network names from prediction files
detect_networks = function(pred_dir, region) {
    # Get all prediction files in the region directory
    if (!dir.exists(pred_dir)) {
        stop("Region directory not found: ", pred_dir)
    }
    
    # List all prediction files
    files = list.files(pred_dir, pattern = sprintf("MODULE_%s_",region), full.names = FALSE)
    files = files[grepl("predicted",files)]
    
    # Extract network names from filenames
    # Pattern: MODULE_region_testing_genotype_network_predicted.txt.gz
    networks = unique(sapply(files, function(f) {
        parts = strsplit(f, "_")[[1]]
        if (length(parts) > 6) {
            return(paste0(parts[4],"_",parts[[6]]))  # Network name is the 5th part
        }else if(length(parts) > 5 & length(parts) <= 6){
          return(paste0(parts[4],"_",parts[[5]]))
        }else if(length(parts) <= 5){
          return(parts[4])
        }
       
    }))
    
    # Remove NULL values and sort
    networks = sort(networks[!sapply(networks, is.null)])
    networks[grepl("dlpfc|caudate|hippo",networks)] = gsub("_","__", networks[grepl("dlpfc|caudate|hippo",networks)])
    
    if (length(networks) == 0) {
        stop("No network names detected in prediction files for region: ", region)
    }
    
    message("Detected networks for region ", region, ": ", paste(networks, collapse = ", "))
    return(networks)
}

# Function to safely read prediction files
safe_read_predictions = function(file_path, model_summaries) {
   library(data.table)
   tmp = fread(file_path)
   tmp = as.data.frame(tmp)
   ## Filter with only training surviving genes
   model_summaries = model_summaries %>% filter(adj_rsq_gene >= 0.01 & pval_gene < 0.05)
   tmp = tmp[,which(colnames(tmp) %in% c("FID","IID",model_summaries$gene))]
   
   return(tmp)
  }

# Set working directory to project root
# setwd("path/to/project/root")
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
    "dplyr", "purrr", "limma","parallel"
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
        # Brain regions to analyze
        regions = c('dlpfc'),
        
        # Directory containing MODULE predictions
        predictions_dir = './predictions/MODULE/',
        
        # List of networks to process
        networks = list.files('./predictions/MODULE/ '),
        
        ## Training model summaries
        model_summary_file = './STEP3.Training_MODULE/output/database/all_networks_summaries',
        
        n_cores            = 10
    ),
    
    output = list(
        testing_genotype = "",
        # Output directory for averaged predictions
        base_dir = './predictions/MODULE/MODULE_averaged/',
        file_prefix = "MODULE_averaged_"
    )
)

# Check environment before proceeding
check_environment(CONFIG)

# Detect networks for each region
CONFIG$input$networks = list()
for (region in CONFIG$input$regions) {
    CONFIG$input$networks[[region]] = detect_networks(CONFIG$input$predictions_dir,region)
}

if(!file.exists(CONFIG$output$base_dir)){dir.create(CONFIG$output$base_dir)}

#' Load and Process Testing Data
#' @param data_path Path to the testing RNA seq data
#' @param tissues Vector of tissue names to analyze
#' @return List of expression data for each tissue
load_testing_data = function(data_path = './rnaseq/your_testing_rnaseq.RData', 
                             tissues = c("dlpfc")) {
    testing_data = get(load(data_path))
    map(tissues %>% set_names(.,.), ~ {
        tissue = .x
        tmp=testing_data[[tissue]]$assays$expression
        rownames(tmp) = paste0(rownames(tmp),"_",rownames(tmp))
        return(tmp)
    })
}

#################################################
## Main Pipeline Function
#################################################

main = function() {
    start_time = Sys.time()
    success_count = 0
    error_count   = 0
    
    # Process each brain region
    for(region in CONFIG$input$regions) {
        message("\nProcessing region: ", region)
        
        tryCatch({
            # Get prediction files for current region
            pred_dir = CONFIG$input$predictions_dir
            if (!dir.exists(pred_dir)) {
                warning("Directory not found: ", pred_dir)
                error_count = error_count + 1
                next
            }
            
            # Get prediction files for each network
            pred_files = list()
            for(network in CONFIG$input$networks[[region]]) {
                pattern = paste0("MODULE_", region)
                files = list.files(pred_dir, pattern = pattern, full.names = TRUE)
                pattern = paste0(network, "_predicted.txt.gz")
                files = list.files(pred_dir, pattern = pattern, full.names = TRUE)
                
                if(length(files) > 0) {
                    pred_files[[network]] = files[1]
                }
            }
            
            if(length(pred_files) == 0) {
                warning("No prediction files found for ", region)
                error_count = error_count + 1
                next
            }
            
            message("Found ", length(pred_files), " prediction files")
            
            ## Additional check: load training summary stats and filter genes  ### TO CHANGE
            # if(CONFIG$output$testing_genotype == "LIBD" && region %in% c("sACC","dACC")){
            #   model_summaries      = readRDS(sprintf("%s_acc.rds", CONFIG$input$model_summary_file))
            # }else{  model_summaries = readRDS(sprintf("%s_%s.rds", CONFIG$input$model_summary_file, region))}
            model_summaries = get(load(sprintf("../../PredictDB/module/MODULE_ExtraDB_lst_%s.RData",region)))
            
            # Read all prediction files
            predictions_list = lapply(pred_files, function(file) {
              safe_read_predictions(file, model_summaries)
            })
            
            predictions_list = predictions_list[!sapply(predictions_list, is.null)]
            
            if(length(predictions_list) == 0) {
                warning("No valid prediction data found for ", region)
                error_count = error_count + 1
                next
            }
            
            message("Successfully loaded ", length(predictions_list), " prediction files")
            
            # Load testing data
            testing_data_true_expr = load_testing_data()
            
            # Get network names
            net.names = names(pred_files)
            
            
            # Average predictions across networks
            averaged_predictions = ensembleAcrossNet(
                model        = "MODULE",
                brain.region = region,
                net.names    = net.names,
                true.expr    = testing_data_true_expr[[region]], 
                outputdir    = CONFIG$output$base_dir,
                outname      = "predicted",
                testing_geno = CONFIG$output$testing_genotype,
                network_predictions = predictions_list,
                pipeline     = "MODULE",
                genestokeep  = NULL,
                n_cores      = CONFIG$input$n_cores
            )
            
            if (is.null(averaged_predictions)) {
                warning("Failed to generate averaged predictions for ", region)
                error_count = error_count + 1
                next
            }
            
            success_count = success_count + 1
            message("Successfully processed ", region)
            
        }, error = function(e) {
            error_count = error_count + 1
            warning("Error processing ", region, ": ", e$message)
        })
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
main()

