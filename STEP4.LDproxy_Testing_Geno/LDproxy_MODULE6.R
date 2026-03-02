# Load required libraries
library(DBI)
library(RSQLite)
library(purrr)
library(dplyr)
library(limma)
library(parallel)
library(tidyr)  
library(magrittr)  
library(data.table)  


# Set working directory to project root
# setwd("path/to/project/root")


###############
## PARAMS
###############

models          = c('MODULE')
regions         = c('dlpfc')
target_file_path = file.path(".", "genotype")  
name_site       = "your_geno_name"
name_site2      = "your_training_dataset_name"
expected_db_count = 48  # Define expected count as a parameter



#########################
## MAIN
#######################


# Iterate over models and regions
walk(models, function(model) {
   
  
  walk(regions, function(brain.region) {
    # Process for the current name_site
    print(paste0(model, " in ", name_site, "; ", toupper(brain.region)))
    
    ## get extras
    extra = readRDS(file.path(".", "STEP3.Training_MODULE", "output", "database", sprintf("all_networks_summaries_%s.rds", brain.region)))
    
    # check if name_site has already expected number of db files
    full_path = file.path(".", "PredictDB_LD", model, name_site, brain.region)
    # Construct the shell command for name_site
    command = paste0("ls ", full_path, " | grep '", name_site, "*.db' | wc -l")
  
    # Execute the shell command and capture the output
    output = system(command, intern = TRUE)
  
    # Convert the output to a numeric value
    file_count = as.numeric(output)
  
    # Print the file count for this name_site
    cat("Checking for:", name_site, "found", file_count, "files.\n")
  
    # Check if the file count matches expected count
    if (file_count != expected_db_count) {
      print(paste0("Creating DB files. Expected ", expected_db_count, " files, found ", file_count))
    } else {
      return(print("already DONE"))
    }
  
  ## Get weights specific to name_site
  weights_file = file.path(".", "PredictDB_LD", model, paste0(model, ".", brain.region, ".extras_weights_", name_site, ".RData"))
  if (!file.exists(weights_file)) {
    stop(paste("Weights file not found:", weights_file))
  }
  weights = get(load(weights_file))
  
  weights         = weights %>% distinct(gene, rsid, network, .keep_all = TRUE) 
  weights$network = gsub(paste0(brain.region,"_network_"), "", weights$network)
  extra$network   = gsub(paste0(brain.region,"_network_"), "", extra$network)
  
  saved_db_file_path = file.path('.', 'PredictDB_LD', model, name_site, brain.region)
  if (!dir.exists(saved_db_file_path)) {
    dir.create(saved_db_file_path, recursive = TRUE)
  }
  
  dbs = unique(extra$network) 
  dbs = gsub(paste0(brain.region,"_network_"), "", dbs)
  
  mclapply(dbs, function(db) {
    # Create DB
    if (!is.na(db)) {
      output_drv_name = paste0("MODULE_", db, '.LD.', name_site2, '.db')
      db_path         = file.path(saved_db_file_path, output_drv_name)
      
      if (file.exists(db_path)) {
        return(print("already DONE"))
      }
      
      tryCatch({
        ## Create DB with modified structure
        con = dbConnect(RSQLite::SQLite(), db_path)
        
        # Create tables with error handling
        dbExecute(con, "
          CREATE TABLE extra (
            extra_id INTEGER PRIMARY KEY,
            description TEXT NOT NULL,
            details TEXT
          );
        ")
        
        dbExecute(con, "
          CREATE TABLE weights (
            weight_id INTEGER PRIMARY KEY,
            value REAL NOT NULL,
            unit TEXT NOT NULL,
            extra_id INTEGER,
            FOREIGN KEY(extra_id) REFERENCES extra(extra_id)
          );
        ")
        
        ## Subset weights and extra for the current db
        weights_net = weights %>% filter(network %in% db)
        extra_net = extra %>% filter(network %in% db, gene_id %in% weights$gene) %>% 
          dplyr::rename(gene = "gene_id", genename = "gene_name", n.snps.in.model = "n_snps_in_model") %>% mutate(pred.perf.R2 = NA, pred.perf.pval=NA, pred.perf.qval=NA)
        
        ## Convert to data.frame if data.table
        if (inherits(weights_net, "data.table")) {
          weights_net = as.data.frame(weights_net)
        }
        if (inherits(extra_net, "data.table")) {
          extra_net = as.data.frame(extra_net)
        }
        
        ## Write tables to the DB
        dbWriteTable(con, 'weights', weights_net, overwrite = TRUE)
        dbWriteTable(con, 'extra', extra_net, overwrite = TRUE)
        
        ## Create indexes
        dbExecute(con, "CREATE INDEX weights_rsid ON weights (rsid)")
        dbExecute(con, "CREATE INDEX weights_gene ON weights (gene)")
        dbExecute(con, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
        dbExecute(con, "CREATE INDEX gene_model_summary ON extra (gene)")
        
        dbDisconnect(con)
        print(paste0("Successfully created database: ", output_drv_name))
      }, error = function(e) {
        print(paste0("Error creating database ", output_drv_name, ": ", e$message))
        if (exists("con") && dbIsValid(con)) {
          dbDisconnect(con)
        }
      })
    }
  }, mc.cores = 1, mc.preschedule = FALSE)
  
  print(paste0(name_site, " DONE"))
  rm(extra)
  rm(weights)
  gc()
  }) 
  gc() # Close walk(regions, function(brain.region)
}) # Close walk(models, function(model)

