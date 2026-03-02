#!/usr/bin/env Rscript

# Load required packages
require(furrr)
require(data.table)
require(purrr)

# Load configuration and utilities
source("config.R")
source("utils.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 1.coTWAS_preprocessing.R <site_index>")
}
site_index <- as.numeric(args[1])

# Set up parallel processing
options(future.globals.maxSize = 20148 * 1024^2)
future::plan("multicore")

# Load metadata and gene sets
metadata_info <- load_metadata(site_index)
site <- metadata_info$site
gene_sets <- load_gene_sets()

# Process each tissue
furrr::future_walk(TISSUES, function(tissue) {
  message("Processing tissue: ", tissue)
  
  # Get predictions
  predictions <- get_predictions(site, tissue)
  
  # Process common genes
  common_predictions <- list(
    E = predictions$epixcan[, colnames(predictions$epixcan) %in% c("IID", gene_sets$common_genes[[tissue]])],
    C = predictions$cis[, colnames(predictions$cis) %in% c("IID", gene_sets$common_genes[[tissue]])],
    M = predictions$module[, colnames(predictions$module) %in% c("IID", gene_sets$common_genes[[tissue]])],
    I = predictions$ingene[, colnames(predictions$ingene) %in% c("IID", gene_sets$common_genes[[tissue]])]
  )
  save_predictions(common_predictions, site, tissue, "common")
  
  # Process unique genes
  unique_predictions <- list(
    E = predictions$epixcan[, colnames(predictions$epixcan) %in% c("IID", gene_sets$unique_genes[[tissue]])],
    C = predictions$cis[, colnames(predictions$cis) %in% c("IID", gene_sets$unique_genes[[tissue]])],
    M = predictions$module[, colnames(predictions$module) %in% c("IID", gene_sets$unique_genes[[tissue]])],
    I = predictions$ingene[, colnames(predictions$ingene) %in% c("IID", gene_sets$unique_genes[[tissue]])]
  )
  save_predictions(unique_predictions, site, tissue, "unique")
})

# Reset parallel processing
future::plan("sequential")

message("Preprocessing completed successfully for site: ", site)
 
