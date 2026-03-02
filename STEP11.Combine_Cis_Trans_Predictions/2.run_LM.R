# Script: 2.run_LM.R
# Purpose: Fits linear models to compare different combinations of predictors (MODULE, INGENE, EpiXcan, CIS)
#          for gene expression prediction in testing data. This script:
#   1. Loads preprocessed testing data predictions
#   2. For each brain region and gene, fits various linear model combinations
#   3. Saves the fitted models for further analysis


###########################
# Load Configuration     #
###########################
# Set working directory to project root
# setwd("path/to/project/root")

# Load configuration
source("./STEP11.Combine_cis_trans_predictions/config.R")

# Set working directory
#setwd(BASE_DIR)

###########################
# Load Required Libraries #
###########################

# Load custom functions
source("./STEP11.Combine_cis_trans_predictions/LM_models_helper_functions.R")

# Data manipulation and statistical analysis libraries
library(limma)  # For strsplit2()
library(data.table)  # For fread()
library(purrr)  # For map functions
library(magrittr)  # For pipe operator %>%

###########################
# Load Preprocessed Data  #
###########################

# Load testing data predictions list from previous step
pred_lst = get(load(file.path(DATA_PATHS$predictions$combined, 
                              FILE_PATTERNS$output$predictions_lst)))

# Extract available brain regions and prediction combinations
regions = BRAIN_REGIONS
combos = names(pred_lst[[1]]$common_predictions)

###########################
# Main Analysis Pipeline  #
###########################

#' Process a single brain region and prediction combination
#' @param region Brain region to process
#' @param combo Prediction combination to analyze
#' @param pred_lst List of predictions
#' @return Fitted models for the region and combination
process_region_combo = function(region, combo, pred_lst) {
  print(paste0("Processing ", region, ": ", combo))
  
  # Get prediction data for current combination
  All.df = data.frame(pred_lst[[region]]$common_predictions[[combo]])
  genes = unique(strsplit2(colnames(All.df), split = '\\.')[, 2])
  
  # Get true expression data
  true.expr_df = pred_lst[[region]][["real.expr"]]
  
  # Fit models for each gene
  sapply(genes, function(gene) {
    print(paste0("Processing ", region, ": ", combo," gene: ",which(genes==gene),"/",length(genes)))
    true.expr = true.expr_df[which(rownames(true.expr_df) %in% rownames(All.df)),]
    true.expr = true.expr[match(rownames(All.df),rownames(true.expr)),]
    stopifnot(identical(rownames(All.df),rownames(true.expr)))
    
    # Prepare data for modeling
    true.expr = true.expr[, gene]
    df.lm = data.frame(cbind(All.df[, grep(gene, colnames(All.df))], true.expr))
    
    # Get column indices for each predictor type
    gene_cols = list(
      MODULE.gene  = grep("MODULE", colnames(df.lm), value = TRUE),
      INGENE.gene  = grep("INGENE", colnames(df.lm), value = TRUE),
      EpiXcan.gene = grep("EpiXcan", colnames(df.lm), value = TRUE),
      CIS.gene     = grep("CIS", colnames(df.lm), value = TRUE)
    )
    
    # Fit appropriate models based on combination type
    fit_linear_models(df.lm, combo, gene_cols)
  }, simplify = FALSE, USE.NAMES = TRUE)
}

# Process each brain region
Common.Predictions.LM = sapply(regions, function(region) {
  if(region=="amygdala"){combos = combos[!grepl("E",combos)]}
  # Process each prediction combination
  model = sapply(combos, function(combo) {
    process_region_combo(region, combo, pred_lst)
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # Save results for current brain region
  save(model, 
       file = file.path(DATA_PATHS$predictions$combined,
                        get_file_path(FILE_PATTERNS$output$model_results,
                                      brain_region = region)))
  
  return(model)
}, simplify = FALSE, USE.NAMES = TRUE)

print("Model fitting completed successfully")


