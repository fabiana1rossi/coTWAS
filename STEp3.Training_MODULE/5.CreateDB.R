#####################################################################
# Create Prediction Database Pipeline
#####################################################################
# This script creates SQLite databases from the training results
# It performs the following steps:
# 1. Processes summary and weights files for each network
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
# setwd("path/to/project/root")

# Source helper functions
source("./STEP3.Training_MODULE/database_helper_functions.R")

#################################################
## Configuration
#################################################
base_dir = getwd()

CONFIG = list(
    summary_dir          = file.path(base_dir, "STEP3.Training_MODULE/output/summary"),
    weights_dir         = file.path(base_dir, "STEP3.Training_MODULE/output/weights"),
    db_dir              = file.path(base_dir, "STEP3.Training_MODULE/output/database"),
    brain_regions       = c('dlpfc'),
    dataset_training    = 'training_dataset_name',
    n_cores             = 1
)

# Create output directories if they don't exist
if (!dir.exists(CONFIG$db_dir)) dir.create(CONFIG$db_dir, recursive = TRUE)

#################################################
## Run Database Creation Pipeline
#################################################

# Process each brain region
for(brain_region in CONFIG$brain_regions) {
    message("Processing brain region: ", brain_region)
    
    # Get network names from summary directory
    network_files = list.files(file.path(CONFIG$summary_dir, brain_region))
    net_names     = gsub("\\.extra\\.txt$", "", network_files)
    
    if(length(net_names) == 0) {
        warning("No network files found for brain region: ", brain_region)
        next
    }
    
    # Create databases for each network
    db_files = makeModelDB1(
        brain.region       = brain_region,
        file_input_summ    = CONFIG$summary_dir,
        file_input_weights = CONFIG$weights_dir,
        drv_file_path      = CONFIG$db_dir,
        net.names          = net_names,
        no_cores           = CONFIG$n_cores,
        dataset_training   = CONFIG$dataset_training
    )
    
    message("Created database files for brain region ", brain_region, ":")
    print(db_files)
}

message("Database creation pipeline completed successfully!") 

#################################################
## Save Combined Summaries and Weights
#################################################

# Initialize lists to store combined data
all_summaries = list()
all_weights   = list()

# Process each brain region
for(brain_region in CONFIG$brain_regions) {
    message("Combining data for brain region: ", brain_region)
    
    # Get network names from summary directory
    network_files = list.files(file.path(CONFIG$summary_dir, brain_region))
    net_names     = gsub("\\.extra\\.txt$", "", network_files)
    
    if(length(net_names) == 0) {
        warning("No network files found for brain region: ", brain_region)
        next
    }
    
    # Create summaries and weights
    summaries_weights = map(net_names, function(net) {
        summary_file = file.path(CONFIG$summary_dir, brain_region, paste0(net, ".extra.txt"))
        weight_file = file.path(CONFIG$weights_dir, brain_region, paste0(net, ".weights.txt"))
        if(file.exists(summary_file) & file.exists(weight_file)) {
            df_s = fread(summary_file)
            ## filter
            df_s = df_s %>% dplyr::filter(adj_rsq_gene_all_data >=0.01 & pval_gene_all_data < .05)
            df_s$network = net
            df_s$brain_region = brain_region
            
            ## weights
            df_w = fread(weight_file)
            df_w = df_w %>% filter(gene %in% df_s$gene_id) %>% mutate(brain_region=brain_region, network=net)
            
            return(list(summary=df_s, weight=df_w))
        }
        return(NULL)
    })
  
    
    # Remove NULL entries and combine
    summaries = do.call(rbind, map(summaries_weights, ~.x[["summary"]]))
    weights   = do.call(rbind, map(summaries_weights, ~.x[["weight"]]))
    
    # Save RData files
    saveRDS(summaries, file.path(CONFIG$db_dir, paste0("all_networks_summaries_",brain_region,".rds")))
    saveRDS(weights, file.path(CONFIG$db_dir, paste0("all_networks_weights_",brain_region,".rds")))
    
    
}

# # Combine all brain regions
# combined_summaries = do.call(rbind, all_summaries)
# combined_weights = do.call(rbind, all_weights)
# 
# # Save RData files
# saveRDS(combined_summaries, file.path(CONFIG$db_dir, "all_networks_summaries.rds"))
# saveRDS(combined_weights, file.path(CONFIG$db_dir, "all_networks_weights.rds"))

message("Successfully saved combined summaries and weights as RData files!") 
