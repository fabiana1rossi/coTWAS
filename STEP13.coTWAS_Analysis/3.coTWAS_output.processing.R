#!/usr/bin/env Rscript

# Load required packages
require(metafor)
require(metap)
require(dplyr)
require(purrr)

# Load configuration and utilities
source("config.R")
source("utils.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 3.coTWAS_output.processing.R <tissue>")
}
tissue <- args[1]

#' Perform meta-analysis for a single gene
#' @param gene_results List of results for a single gene across sites
#' @return List of meta-analysis results
perform_meta_analysis <- function(gene_results) {
  # Filter valid results
  valid_results <- Filter(length, gene_results)
  
  if (length(valid_results) > 1) {
    sites <- names(valid_results)
    predictors <- sapply(valid_results, function(site) site$predictors)
    
    # Meta-analysis for linear regression
    eff_size <- as.numeric(sapply(valid_results, function(site) site$gene_df_linear$statistic[2]))
    n_size <- sapply(valid_results, function(site) site$df_stats_linear$nobs[2])
    n_predictors <- as.numeric(sapply(valid_results, function(site) nrow(site$gene_df_linear) - 1))
    
    meta_input <- escalc(measure = "ZPCOR", ti = eff_size, ni = n_size, mi = n_predictors)
    meta_linear <- rma.uni(yi = meta_input$yi,
                          vi = meta_input$vi,
                          control = list(stepadj = 0.5, maxiter = 10000),
                          method = "REML")
    
    # Meta-analysis for logistic regression
    eff_size <- as.numeric(sapply(valid_results, function(site) site$gene_df_logistic$estimate[2]))
    se <- as.numeric(sapply(valid_results, function(site) site$gene_df_logistic$std.error[2]))
    
    meta_logistic <- rma.uni(yi = as.numeric(eff_size),
                            sei = as.numeric(se),
                            control = list(stepadj = 0.5, maxiter = 10000),
                            method = "REML")
    
    # Meta-analysis for p-values
    p_values <- as.numeric(sapply(valid_results, function(site) site$anova_p_logistic))
    p_values[is.na(p_values)] <- 1
    meta_p_stouffer <- as.numeric(metap::sumz(p_values, weights = sqrt(n_size))$p)
    
    # Calculate pseudo R-squared
    pseudoR <- as.numeric(sapply(valid_results, function(site) site$pseudoR))
    
    return(list(
      sites = sites,
      predictors = predictors,
      meta_linear_beta = as.numeric(meta_linear$beta),
      meta_linear_se = as.numeric(meta_linear$se),
      meta_linear_p = meta_linear$pval,
      meta_linear_QEp = as.numeric(meta_linear$QEp),
      meta_logistic_beta = as.numeric(meta_logistic$beta),
      meta_logistic_se = as.numeric(meta_logistic$se),
      meta_logistic_p = meta_logistic$pval,
      meta_logistic_QEp = as.numeric(meta_logistic$QEp),
      pred_linear_p = as.numeric(sapply(valid_results, function(site) site$gene_df_linear$p.value[2])),
      pred_logistic_p = as.numeric(sapply(valid_results, function(site) site$gene_df_logistic$p.value[2])),
      pseudoR = pseudoR,
      meta_p_stouffer = meta_p_stouffer,
      n_size = n_size
    ))
  } else if (length(valid_results) == 1) {
    sites <- names(valid_results)
    predictors <- sapply(valid_results, function(site) site$predictors)
    
    # Extract single site results
    site <- valid_results[[1]]
    return(list(
      sites = sites,
      predictors = predictors,
      meta_linear_beta = NA,
      meta_linear_se = NA,
      meta_linear_p = NA,
      meta_linear_QEp = NA,
      meta_logistic_beta = NA,
      meta_logistic_se = NA,
      meta_logistic_p = NA,
      meta_logistic_QEp = NA,
      pred_linear_p = as.numeric(site$gene_df_linear$p.value[2]),
      pred_logistic_p = as.numeric(site$gene_df_logistic$p.value[2]),
      pseudoR = site$pseudoR,
      meta_p_stouffer = NA,
      n_size = site$df_stats_linear$nobs[2]
    ))
  } else {
    return(list(
      sites = NA,
      predictors = NA,
      meta_linear_beta = NA,
      meta_linear_se = NA,
      meta_linear_p = NA,
      meta_linear_QEp = NA,
      meta_logistic_beta = NA,
      meta_logistic_se = NA,
      meta_logistic_p = NA,
      meta_logistic_QEp = NA,
      pred_linear_p = NA,
      pred_logistic_p = NA,
      pseudoR = NA,
      meta_p_stouffer = NA,
      n_size = NA
    ))
  }
}

# Load results for the tissue
results_files <- list.files(OUTPUT_PATH)
results_tissue <- grep(tissue, results_files, value = TRUE)

# Load and process results
results <- purrr::map(results_tissue, function(x) {
  load(file.path(OUTPUT_PATH, x))
  return(results_coTWAS)
})
names(results) <- limma::strsplit2(results_tissue, "_")[, 1]

# Transpose results to get genes as first level
results <- purrr::transpose(results)

# Get information about available sites
check <- sapply(results, function(gene) sum(sapply(gene, function(site) length(site) == 11)))
check_miss <- sapply(results, function(gene) names(Filter(function(x) length(x) != 11, gene)), simplify = FALSE)
check_in <- sapply(results, function(gene) names(Filter(function(x) length(x) == 11, gene)), simplify = FALSE)

# Perform meta-analysis for each gene
meta_results <- purrr::imap(results, function(gene, name) {
  perform_meta_analysis(gene)
})

# Create output data frame
output_coTWAS <- data.frame(
  tissue = tissue,
  gene = names(results),
  N_site = as.numeric(check),
  sites = I(sapply(meta_results, function(x) x$sites)),
  pvalues = I(sapply(meta_results, function(x) x$pred_logistic_p)),
  predictors = I(sapply(meta_results, function(x) x$predictors)),
  pseudoR2_average = as.numeric(sapply(meta_results, function(x) mean(x$pseudoR))),
  pseudoR2_median = as.numeric(sapply(meta_results, function(x) median(x$pseudoR))),
  beta_linear = as.numeric(sapply(meta_results, function(x) x$meta_linear_beta)),
  se_linear = as.numeric(sapply(meta_results, function(x) x$meta_linear_se)),
  beta_logistic = as.numeric(sapply(meta_results, function(x) x$meta_logistic_beta)),
  se_logistic = as.numeric(sapply(meta_results, function(x) x$meta_logistic_se)),
  meta_linear_p = as.numeric(sapply(meta_results, function(x) x$meta_linear_p)),
  meta_logistic_p = as.numeric(sapply(meta_results, function(x) x$meta_logistic_p)),
  meta_linear_QEp = as.numeric(sapply(meta_results, function(x) x$meta_linear_QEp)),
  meta_logistic_QEp = as.numeric(sapply(meta_results, function(x) x$meta_logistic_QEp))
)

# Add gene information
output_coTWAS <- add_gene_info(output_coTWAS)

# Save results
save(output_coTWAS,
     file = file.path(OUTPUT_PATH, paste0("output_", tissue, ".RData")))

message("Processing completed successfully for tissue: ", tissue)
