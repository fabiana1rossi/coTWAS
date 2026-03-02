# Utility functions for coTWAS analysis

#' Load and validate metadata
#' @param site_index Index of the site to analyze
#' @return List containing metadata and site name
load_metadata <- function(site_index) {
  if (!file.exists(METADATA_PATH)) {
    stop("Metadata file not found at: ", METADATA_PATH)
  }
  
  load(METADATA_PATH)
  site <- names(metadata.PRS.sub)[site_index]
  
  if (is.null(site)) {
    stop("Invalid site index: ", site_index)
  }
  
  return(list(
    metadata = metadata.PRS.sub,
    site = site
  ))
}

#' Load gene sets
#' @return List containing common and unique genes
load_gene_sets <- function() {
  if (!file.exists(COMMON_GENES_PATH) || !file.exists(UNIQUE_GENES_PATH)) {
    stop("Gene set files not found")
  }
  
  load(COMMON_GENES_PATH)
  load(UNIQUE_GENES_PATH)
  
  return(list(
    common_genes = Common_Genes$Chosen.genes,
    unique_genes = Unique_Genes$Chosen.genes
  ))
}

#' Get prediction files for a specific site and tissue
#' @param site Site name
#' @param tissue Tissue name
#' @return List of prediction data frames
get_predictions <- function(site, tissue) {
  # Get file patterns
  epixcan_files <- grep("predicted.txt", list.files(EPIXCAN_PATH), value = TRUE)
  cis_files <- grep("predicted.txt", list.files(CIS_PATH), value = TRUE)
  module_files <- grep("MODULE", list.files(MODULE_PATH), value = TRUE)
  ingene_files <- grep("INGENE", list.files(INGENE_PATH), value = TRUE)
  
  # Load predictions
  df_epixcan <- fread(file.path(EPIXCAN_PATH, grep(site, epixcan_files, value = TRUE)))
  df_cis <- fread(file.path(CIS_PATH, grep(paste0(tissue, "_", site), cis_files, value = TRUE)))
  df_module <- fread(file.path(MODULE_PATH, grep(paste0(site, "_", tissue), module_files, value = TRUE)))
  df_ingene <- fread(file.path(INGENE_PATH, grep(paste0(site, "_", tissue), ingene_files, value = TRUE)))
  
  # Standardize column names
  colnames(df_module)[1] <- colnames(df_ingene)[1] <- "IID"
  
  return(list(
    epixcan = df_epixcan,
    cis = df_cis,
    module = df_module,
    ingene = df_ingene
  ))
}

#' Save predictions to RData file
#' @param predictions List of prediction data frames
#' @param site Site name
#' @param tissue Tissue name
#' @param gene_set_type Type of gene set ("common" or "unique")
save_predictions <- function(predictions, site, tissue, gene_set_type) {
  output_file <- file.path(INPUT_PATH, paste0(site, "_", tissue, "_predictions.", gene_set_type, ".RData"))
  save(predictions, file = output_file)
}

#' Check if a gene is in the MHC region
#' @param seqnames Chromosome name
#' @param start Start position
#' @param end End position
#' @return Logical indicating if gene is in MHC region
is_in_mhc <- function(seqnames, start, end) {
  seqnames == MHC_REGION$chr & 
    start >= MHC_REGION$start & 
    end <= MHC_REGION$end
}

#' Add gene information to results
#' @param results Data frame of results
#' @return Data frame with added gene information
add_gene_info <- function(results) {
  if (!file.exists(GENE_INFO_PATH)) {
    stop("Gene info file not found at: ", GENE_INFO_PATH)
  }
  
  gene_info <- readRDS(GENE_INFO_PATH)
  
  results$genes.symbol <- gene_info$Symbol.x[match(results$gene, gene_info$ensemblID)]
  results$seqnames <- gene_info$seqnames[match(results$gene, gene_info$ensemblID)]
  results$start <- gene_info$start[match(results$gene, gene_info$ensemblID)]
  results$end <- gene_info$end[match(results$gene, gene_info$ensemblID)]
  results$MHC <- ifelse(is_in_mhc(results$seqnames, results$start, results$end), "MHC", "noMHC")
  results$gene_type <- gene_info$gene_type[match(results$gene, gene_info$ensemblID)]
  
  return(results)
} 