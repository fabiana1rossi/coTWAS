# This script processes and combines gene expression predictions from two different INGENE models:
# 1. INGENE_EpiXcan: A model that uses epigenetic data for prediction
# 2. INGENE_CIS: A model that uses cis-regulatory information for prediction
# The script selects the best performing model for each gene and combines the predictions
# This script uses the averaged network predictions from 1.Average_net_INGENE.R and
# the selected best models from 2.Select_best_cis_model.R

# Set working directory and base directory
# setwd("path/to/project/root")


# Load required libraries
require(purrr)      # For functional programming
library(tidyverse)  # For data manipulation and visualization
library(IRanges)    # For genomic range operations
library(limma)      # For microarray data analysis
library(data.table)

############
###### PARAMS
#######
# Directory containing INGENE results for GTEx data
result_dir = "./predictions/INGENE"
# Directory containing averaged network predictions
avg_pred_dir = sprintf("%s/INGENE_averaged", result_dir)
# Directory containing selected model performance
model_perf_dir = "./predictions/"

###############
## FUNCTIONS
###############
#' Load and Process Testing Data
#' @param data_path Path to the testing RNA seq data
#' @param tissues Vector of tissue names to analyze
#' @return List of expression data for each tissue
load_testing_data = function(data_path = './rnaseq/your_rnaseq.RData', 
                             tissues = c("dlpfc")) {
  testing_data = get(load(data_path))
  map(tissues %>% set_names(.,.), ~ {
    tissue = .x
    if(tissue %in% c("dACC","sACC")){tissue="acc"}
    tmp=testing_data[[tissue]][['assays']][['expression']]
    rownames(tmp) = paste0(rownames(tmp),"_",rownames(tmp))
    tmp
  })
}

# Function to process and combine imputed gene expression data from both INGENE models
# Parameters:
# - region: Brain region being analyzed
# - ingene_epixcan_gtex.imputed.data: Imputed expression data from EpiXcan model
# - ingene_epi_perf: Performance metrics for EpiXcan model
# - ingene_cis_gtex.imputed.data: Imputed expression data from CIS model
# - ingene_cis_perf: Performance metrics for CIS model
# - cov_df: Optional covariates data frame
getImputedData = function(real_df, region, ingene_epixcan_gtex.imputed.data, ingene_epi_perf, 
                         ingene_cis_gtex.imputed.data, ingene_cis_perf, cov_df = NULL){
  
  # Filter performance data to only include Averaged Ensemble network results
  ingene_epi_perf       = ingene_epi_perf[ingene_epi_perf$network=="Averaged_Ensemble",]
  ingene_epi_perf$Model = "INGENE_EpiXcan"
  
  ingene_cis_perf       = ingene_cis_perf[ingene_cis_perf$network=="Averaged_Ensemble",]
  ingene_cis_perf$Model = "INGENE_CIS"
  
  # Compare performance between models and select the best one for each gene
  # Join performance metrics from both models and determine which performs better
  epi.cis_perf  = inner_join(ingene_epi_perf, ingene_cis_perf, by = c("gene")) %>% 
    dplyr::rename(r_square_adjusted_epi = "r_square_adjusted.x", r_square_adjusted_cis = "r_square_adjusted.y") %>%
    mutate(win_model = case_when(r_square_adjusted_epi > r_square_adjusted_cis ~ "EpiXcan",
                                 r_square_adjusted_cis > r_square_adjusted_epi ~ "CIS",
                                 r_square_adjusted_epi == r_square_adjusted_cis ~ "Same"))
  
  # Print statistics about model selection
  print(length(epi.cis_perf$gene))
  print(table(epi.cis_perf$win_model))
  
  # Filter performance dataframes to keep only the winning model for each gene
  ingene_cis_perf  = ingene_cis_perf %>% filter(!gene %in% epi.cis_perf$gene[epi.cis_perf$win_model=="EpiXcan"])
  ## IF Same choose EpiXcan
  ingene_cis_perf  = ingene_cis_perf %>% filter(!gene %in% epi.cis_perf$gene[epi.cis_perf$win_model=="Same"])
  ingene_epi_perf   = ingene_epi_perf %>% filter(!gene %in% epi.cis_perf$gene[epi.cis_perf$win_model=="CIS"])
  
  # Save combined performance metrics
  save(ingene_cis_perf, file = paste0("./predictions/INGENE/INGENE_averaged/CIS/INGENE_your_testing_genotype_name_performance_",region,"_CIS.RData"))
  save(ingene_epi_perf, file = paste0("./predictions/INGENE/INGENE_averaged/EpiXcan/INGENE_your_testing_genotype_name_performance_",region,"_EpiXcan.RData"))
  
  # Combine performance metrics from both models
  ingene_total_perf = rbind(ingene_cis_perf,ingene_epi_perf)
  
  # Save combined performance metrics
  save(ingene_total_perf, file = paste0("./predictions/INGENE/INGENE_averaged/INGENE_your_testing_genotype_name_performance_",region,"_total.RData"))
  
  # Subset imputed expression data to include only genes from the winning model
  # For EpiXcan model
  ingene_epixcan_gtex.imputed.data   = ingene_epixcan_gtex.imputed.data[,which(colnames(ingene_epixcan_gtex.imputed.data) %in% c("FID","IID",ingene_epi_perf$gene))]
  ingene_epixcan_gtex.imputed.data$FID = ingene_epixcan_gtex.imputed.data$IID = rownames(ingene_epixcan_gtex.imputed.data)
  
  # For CIS model
  ingene_cis_gtex.imputed.data      = ingene_cis_gtex.imputed.data[,which(colnames(ingene_cis_gtex.imputed.data) %in% c("FID","IID",ingene_cis_perf$gene))]
  ingene_cis_gtex.imputed.data$FID  =  ingene_cis_gtex.imputed.data$IID = rownames(ingene_cis_gtex.imputed.data)
  
  # Verify sample matching between models
  stopifnot(identical(rownames(ingene_cis_gtex.imputed.data),rownames(ingene_epixcan_gtex.imputed.data)))
  
  # Combine predictions from both models into a single dataframe
  INGENE_predicted = cbind(ingene_epixcan_gtex.imputed.data,ingene_cis_gtex.imputed.data)
  INGENE_predicted$FID = INGENE_predicted$IID = rownames(INGENE_predicted)
  
  ## Subset with region samples
  INGENE_predicted  = INGENE_predicted[which(INGENE_predicted$IID %in% rownames(real_df)),]
  
  
  # Save combined predictions
  write.table(INGENE_predicted, file = paste0("./predictions/INGENE/INGENE_averaged/",region,"_AveragedNetworks_predicted.txt"))

  return(print("DONE"))
}

# Process data for multiple brain regions
# Regions: dlpfc (dorsolateral prefrontal cortex), dACC (dorsal anterior cingulate cortex),
# sACC (subgenual anterior cingulate cortex), hippo (hippocampus), caudate, and amygdala
walk(c("dACC","sACC","dlpfc") %>% set_names(.,.), function(region){
  
  print(region)
  
  # Load real GTEx expression data for the current region
  gtex_true_expr = load_testing_data()
  #gtex_true_expr$sACC = gtex_true_expr$dACC = gtex_true_expr$acc
  
  gtex.real.data=gtex_true_expr[[region]]
  
  # Load averaged network predictions from CIS INGENE model
  ingene_cis_gtex.imputed.data  = fread(file.path(avg_pred_dir, "CIS",
                                                  paste0(region,"_CIS_AveragedNetworks_predicted.txt")))
  ingene_cis_gtex.imputed.data  = as.data.frame(ingene_cis_gtex.imputed.data)
  rownames(ingene_cis_gtex.imputed.data) = ingene_cis_gtex.imputed.data[[1]]
  ingene_cis_gtex.imputed.data[[1]] = NULL
  
  if(region != "amygdala"){
    # Load averaged network predictions from EpiXcan INGENE model
    ingene_epixcan_gtex.imputed.data  = fread(file.path(avg_pred_dir, "EpiXcan",
                                                        paste0(region,"_EpiXcan_AveragedNetworks_predicted.txt")))
    ingene_epixcan_gtex.imputed.data  = as.data.frame(ingene_epixcan_gtex.imputed.data)
    rownames(ingene_epixcan_gtex.imputed.data) = ingene_epixcan_gtex.imputed.data[[1]]
    ingene_epixcan_gtex.imputed.data[[1]] = NULL
    
    
    
    # Load selected model performance metrics
    ingene_epi_perf = get(load(file.path(avg_pred_dir, "EpiXcan",
                                         paste0("INGENE_your_testing_genotype_name_performance_",region,".RData"))))
    ingene_cis_perf = get(load(file.path(avg_pred_dir, "CIS",
                                         paste0("INGENE_your_testing_genotype_name_performance_",region,".RData"))))
    
    
    # Process and combine predictions from both models
    data = getImputedData(real_df                          = gtex.real.data, 
                          region                           = region,
                          ingene_epixcan_gtex.imputed.data = ingene_epixcan_gtex.imputed.data, 
                          ingene_epi_perf                  = ingene_epi_perf,
                          ingene_cis_gtex.imputed.data     = ingene_cis_gtex.imputed.data,
                          ingene_cis_perf                  = ingene_cis_perf,
                          cov_df                           = NULL)
  }else{
    
    ingene_cis_perf = get(load(file.path(avg_pred_dir, "CIS",
                                         paste0("INGENE_your_testing_genotype_name_performance_",region,".RData"))))
    # Save combined performance metrics
    save(ingene_cis_perf, file = paste0("./predictions/INGENE/INGENE_averaged/INGENE_your_testing_genotype_name_performance_",region,"_total.RData"))
    save(ingene_cis_perf, file = paste0("./predictions/INGENE/INGENE_averaged/CIS/INGENE_your_testing_genotype_name_performance_",region,"_CIS.RData"))
    
    
    ingene_cis_gtex.imputed.data$FID = ingene_cis_gtex.imputed.data$IID = rownames(ingene_cis_gtex.imputed.data)
    
    ## Subset with region samples
    INGENE_predicted  = ingene_cis_gtex.imputed.data[which(ingene_cis_gtex.imputed.data$IID %in% rownames(gtex.real.data)),]
    
    
    # Save combined predictions
    write.table(INGENE_predicted, file = paste0("./predictions/INGENE/INGENE_averaged/",region,"_AveragedNetworks_predicted.txt"))
    
  }


})
ls



