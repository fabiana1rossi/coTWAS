#!/usr/bin/env Rscript

# Load required libraries
library(DBI)
library(RSQLite)
library(purrr)
library(dplyr)
library(limma)
library(parallel)
library(tidyr)  # Added for data manipulation
library(magrittr)  # Added for pipe operations
library(data.table)  # Added for data.table support


# Set working directory to project root
# setwd("path/to/project/root")


###############
## PARAMS
###############

models          = c('CIS')
regions         = c('dlpfc')
target_file_path = file.path(".", "genotype")  
name_site          = "your_genotype_name"     # GTeX.v9.phased.geno.maf.bfile.updated.EUR   CMC_EUR.IBD.geno0.02.hwe10e6.mind0.02.maf0.01.hg20  LIBD_TopMed.geno.maf.EUR.chr1_22.noindels.overlap.GTEX
name_site2         = "suffix_output_file"




#########################
## MAIN
#######################


# Iterate over models and regions
walk(models, function(model) {
  
  walk(regions, function(brain.region) {
    # Process for the current name_site
    print(paste0(model, " in ", name_site, "; ", toupper(brain.region)))
    
    ## Get extras
    extra = read.csv(paste0("./STEP2.Training_cis/output/database/extra_",brain.region,".csv"))

    ## Get weights specific to name_site
    weights_file = file.path(".", "PredictDB_LD", model, paste0(model, ".", brain.region, ".extras_weights_", name_site, ".RData"))
    if (!file.exists(weights_file)) {
      stop(paste("Weights file not found:", weights_file))
    }
    weights = get(load(weights_file))
    
    weights = weights %>% distinct(gene, rsid, .keep_all = TRUE) 
    
    saved_db_file_path = file.path('.', 'PredictDB_LD', model, name_site, brain.region)
    if (!dir.exists(saved_db_file_path)) {
      dir.create(saved_db_file_path, recursive = TRUE)
    }
  
   
    output_drv_name = paste0(model,".LD.", name_site2, '.db')
    db_path         = file.path(saved_db_file_path, output_drv_name)
      
    if (file.exists(db_path)) {
      return(print("already DONE"))
    }
    
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
    extra = extra %>% filter(gene %in% weights$gene) 
    ## Convert to data.frame if data.table
    if (inherits(weights, "data.table")) {
      weights = as.data.frame(weights)
    }
    if (inherits(extra, "data.table")) {
      extra = as.data.frame(extra)
    }
      
    ## Write tables to the DB
    dbWriteTable(con, 'weights', weights, overwrite = TRUE)
    dbWriteTable(con, 'extra', extra, overwrite = TRUE)
      
    ## Create indexes
    dbExecute(con, "CREATE INDEX weights_rsid ON weights (rsid)")
    dbExecute(con, "CREATE INDEX weights_gene ON weights (gene)")
    dbExecute(con, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
    dbExecute(con, "CREATE INDEX gene_model_summary ON extra (gene)")
      
    dbDisconnect(con)
    print(paste0("Successfully created database: ", output_drv_name))
    
    print(paste0(name_site, " DONE"))
    rm(extra)
    rm(weights)
    gc()
  }) 
  gc() # Close walk(regions, function(brain.region)
}) # Close walk(models, function(model)

