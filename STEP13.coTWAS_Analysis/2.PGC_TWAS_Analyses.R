#!/usr/bin/env Rscript

# Load required packages
require(rcompanion)
require(broom)
require(dplyr)
require(mgsub)
require(purrr)

# Load configuration and utilities
source("config.R")
source("utils.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 2.PGC_TWAS_Analyses.R <tissue> <site_index>")
}
tissue <- args[1]
site_index <- as.numeric(args[2])

# Load metadata and gene sets
metadata_info <- load_metadata(site_index)
site <- metadata_info$site
df <- data.frame(metadata_info$metadata[[site]])
df_ncol <- ncol(df)

# Load gene sets
gene_sets <- load_gene_sets()
genes <- purrr::map2(gene_sets$common_genes, gene_sets$unique_genes, function(x, y) rbind(x, y))

# Load predictions
load(file.path(INPUT_PATH, paste0(site, "_", tissue, "_predictions.common.RData")))
predictions_common <- predictions
load(file.path(INPUT_PATH, paste0(site, "_", tissue, "_predictions.unique.RData")))
predictions_unique <- predictions
predictions <- purrr::map2(predictions_common, predictions_unique, function(x, y) merge(x, y, by = "IID"))

# Load model fits
model_file <- file.path(OUTPUT_BASE_PATH, paste0("Common.Predictions.", tissue, ".RData"))
if (!file.exists(model_file)) {
  stop("Model file not found: ", model_file)
}
load(model_file)

# Process model names
model <- lapply(model, purrr::transpose)
names(model) <- mgsub::mgsub(names(model),
                            c("MIE_common_prediction", "MI_common_prediction", "ME_common_prediction",
                              "IE_common_prediction", "MIC_common_prediction", "MC_common_prediction",
                              "IC_common_prediction"),
                            c("EMI", "MI", "EM", "EI", "CMI", "CM", "CI"))

model <- lapply(model, function(x) {
  names(x) <- mgsub::mgsub(names(x),
                          c("MIE", "MI", "ME", "IE", "MIC", "MC", "IC"),
                          c("EMI", "MI", "EM", "EI", "CMI", "CM", "CI"))
  x
})

# Clean up environment
rm(list = setdiff(ls(), c("genes", "df", "df_ncol", "site", "tissue", "predictions", "model")))

#' Run analysis for a single gene
#' @param gene Gene name
#' @param tissue Tissue name
#' @param df Data frame with metadata
#' @param df_ncol Number of columns in original data frame
#' @param predictions List of predictions
#' @param model List of models
#' @return List of analysis results
analyze_gene <- function(gene, tissue, df, df_ncol, predictions, model) {
  # Get model information
  mm_name <- genes[[tissue]]$model[genes[[tissue]]$gene %in% gene]
  models <- limma::strsplit2(mm_name, "")[1,]
  
  # Prepare data
  df[, models] <- sapply(models, function(x) predictions[[x]][[gene]][match(paste0(df$FID, "_", df$IID), predictions[[x]]$IID)], USE.NAMES = FALSE, simplify = FALSE)
  df <- df[, apply(df, 2, function(x) !all(is.na(x)))]
  models_available <- colnames(df)[(df_ncol + 1):ncol(df)]
  check_df <- ncol(df) != df_ncol
  check_m <- paste0(models_available, collapse = "") == mm_name
  check_length_m <- length(models_available) > 1
  
  # Prepare covariates and outcome
  df$Dx <- factor(ifelse(df$Dx == 1, 0, 1))
  df$Sex <- factor(df$Sex)
  cov <- c("Sex", grep("^C[0-9]", colnames(df), value = TRUE))
  outcome <- "Dx"
  
  if (!check_df) {
    return(NULL)
  }
  
  # Prepare predictor based on conditions
  if (check_m && check_length_m) {
    # Predict model from GTEx
    colnames(df) <- mgsub::mgsub(colnames(df),
                                c("^E$", "^M$", "^I$", "^C$"),
                                c(paste0("EpiXcan.", gene),
                                  paste0("MODULE.", gene),
                                  paste0("INGENE.", gene),
                                  paste0("CIS.", gene)))
    mm <- model[[mm_name]][[paste0("model.", mm_name)]][[gene]]
    predictors <- colnames(mm$model)[-1]
    df$predictor <- predict(mm, df[, predictors])
    predictors <- mm_name
  } else if (check_m && !check_length_m) {
    predictors <- mm_name
    df$predictor <- df[[predictors]]
  } else if (!check_m && check_length_m) {
    # Predict model from GTEx
    colnames(df) <- mgsub::mgsub(colnames(df),
                                c("^E$", "^M$", "^I$", "^C$"),
                                c(paste0("EpiXcan.", gene),
                                  paste0("MODULE.", gene),
                                  paste0("INGENE.", gene),
                                  paste0("CIS.", gene)))
    mm_name_available <- paste0(models_available, collapse = "")
    mm <- model[[mm_name]][[paste0("model.", mm_name_available)]][[gene]]
    predictors <- colnames(mm$model)[-1]
    df$predictor <- predict(mm, df[, predictors])
    predictors <- mm_name_available
  } else if (!check_m && !check_length_m) {
    predictors <- paste0(models_available, collapse = "")
    df$predictor <- df[[predictors]]
  }
  
  # Run linear regression
  fit_model <- paste0(c("predictor", paste0(c(outcome, cov), collapse = "+")), collapse = "~")
  gene_fit <- lm(as.formula(fit_model), data = df)
  gene_df_linear <- tidy(gene_fit)
  gene_stats_linear <- glance(gene_fit)
  
  null_model <- paste0(c("predictor", paste0(cov, collapse = "+")), collapse = "~")
  null_fit <- lm(as.formula(null_model), data = df)
  null_df_linear <- tidy(null_fit)
  null_stats_linear <- glance(null_fit)
  
  anova_p_linear <- anova(null_fit, gene_fit)$`Pr(>F)`[2]
  
  # Run logistic regression
  null_model <- paste0(c(outcome, paste0(cov, collapse = "+")), collapse = "~")
  null_fit <- glm(as.formula(null_model), data = df, family = "binomial")
  null_df_logistic <- tidy(null_fit)
  null_stats_logistic <- glance(null_fit)
  
  fit_model <- paste0(c(outcome, paste0(c("predictor", cov), collapse = "+")), collapse = "~")
  gene_fit <- glm(as.formula(fit_model), data = df, family = "binomial")
  gene_df_logistic <- tidy(gene_fit)
  gene_stats_logistic <- glance(gene_fit)
  
  anova_p_logistic <- anova(null_fit, gene_fit, test = "LRT")$`Pr(>Chi)`[2]
  
  # Calculate Nagelkerke R2
  fit_nagelkerke <- nagelkerke(fit = gene_fit, null = null_fit)
  pseudoR <- fit_nagelkerke[[2]][3, 1]
  pseudoR_pval <- fit_nagelkerke[[3]][4]
  
  return(list(
    predictors = predictors,
    gene_df_linear = gene_df_linear,
    gene_df_logistic = gene_df_logistic,
    null_df_linear = null_df_linear,
    null_df_logistic = null_df_logistic,
    df_stats_linear = rbind(null_stats_linear, gene_stats_linear) %>% mutate(models = c("null.fit", "gene.fit")),
    df_stats_logistic = rbind(null_stats_logistic, gene_stats_logistic) %>% mutate(models = c("null.fit", "gene.fit")),
    pseudoR = pseudoR,
    pseudoR_pval = pseudoR_pval,
    anova_p_linear = anova_p_linear,
    anova_p_logistic = anova_p_logistic
  ))
}

# Run analysis for all genes
results_coTWAS <- sapply(genes[[tissue]]$gene,
                        function(gene) analyze_gene(gene, tissue, df, df_ncol, predictions, model),
                        USE.NAMES = TRUE,
                        simplify = FALSE)

# Save results
save(results_coTWAS,
     file = file.path(OUTPUT_PATH, paste0(site, "_", tissue, "_coTWAS.Analysis.RData")))

message("Analysis completed successfully for site: ", site, " and tissue: ", tissue)
