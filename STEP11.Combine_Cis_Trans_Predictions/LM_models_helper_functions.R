# Script: LM_models_helper_functions.R
# Purpose: Contains helper functions for linear modeling analysis of gene expression predictions.
#          This file provides functions for:
#   1. Data preparation and manipulation
#   2. Expression data processing
#   3. Model prediction handling
#   4. Statistical analysis utilities
# Author: Original by Fabiana, enhanced documentation added



###########################
# Data Manipulation Utils #
###########################

#' Find Common Rownames Across Data Frames
#' @description Finds common rownames based on the dataframe with the least number of rows
#' @param df_list List of data frames to compare
#' @return Vector of common rownames
find_common_rownames = function(df_list) {
  # Find the dataframe with the fewest rows
  min_rows_df = df_list[[which.min(map_int(df_list, nrow))]]
  
  # Get the rownames of the dataframe with the fewest rows
  common_rownames = rownames(min_rows_df)
  
  return(common_rownames)
}

#' Match Rownames Across Data Frames
#' @description Subsets all dataframes to include only common rows
#' @param df_list List of data frames to match
#' @return List of data frames with matched rows
match_rownames = function(df_list) {
  common_rownames = find_common_rownames(df_list)
  matched_df_list = map(df_list, ~ .x[common_rownames, ])
  return(matched_df_list)
}

#' Match Both Rownames and Colnames Across Data Frames
#' @description Subsets all dataframes to include only common rows and columns
#' @param df_list List of data frames to match
#' @return List of data frames with matched rows and columns
match_rownames_colnames = function(df_list) {
  common_rownames = find_common_rownames(df_list)
  matched_df_list = map(df_list, ~ .x[common_rownames, ])
  return(matched_df_list)
}

###########################
# Data Processing Utils   #
###########################

#' Process MODULE Data
#' @description Process and prepare MODULE prediction data
#' @param dataset Name of the dataset
#' @param input.dir Input directory path
#' @param brain.region Brain region name
#' @return List containing processed MODULE data and performance metrics
process_module_data = function(dataset, input.dir, brain.region, samples) {
  # Construct file paths
  pred_file = file.path(input.dir, "MODULE/MODULE_averaged",
                       get_method_file_path("module", "predictions",
                                          brain_region = brain.region))
  perf_file = file.path("./predictions/LIBD_training/GTeX.v9/MODULE/MODULE_averaged",
                       get_method_file_path("module", "performance",
                                          dataset = "GTEx",
                                          brain_region = brain.region))
  
  # Read and process data
  predicted_expr = fread(pred_file) %>% as.data.frame()
  rownames(predicted_expr) = predicted_expr$FID
  predicted_expr[[1]] = NULL
  predicted_expr$FID = predicted_expr$IID = NULL
  
  predicted_expr = predicted_expr[which(rownames(predicted_expr) %in% samples),]
  
  # Get performance predictions
  perf_module = get(load(perf_file))
  
  # Subset predictions
  predicted_expr = predicted_expr[, which(colnames(predicted_expr) %in% perf_module$gene[perf_module$network == "Averaged_Ensemble"])]
  
  return(list(
    predictions = predicted_expr,
    performance = perf_module,
    samples = rownames(predicted_expr)
  ))
}

#' Process EpiXcan Data
#' @description Process and prepare EpiXcan prediction data
#' @param dataset Name of the dataset
#' @param input.dir Input directory path
#' @param brain.region Brain region name
#' @return List containing processed EpiXcan data and performance metrics
process_epixcan_data = function(dataset, input.dir, brain.region, samples) {
  # Construct file paths
  pred_file = file.path(input.dir, "EpiXcan",
                       get_method_file_path("epixcan", "predictions",
                                          brain_region = brain.region,
                                          dataset = dataset))
  perf_file = file.path("./predictions/LIBD_training/GTeX.v9/EpiXcan",
                       get_method_file_path("epixcan", "performance",
                                          dataset = "GTEx",
                                          brain_region = brain.region))
  
  # Read and process data
  predicted_expr = fread(pred_file) %>% as.data.frame()
  rownames(predicted_expr) = predicted_expr$FID
  predicted_expr[[1]] = NULL
  predicted_expr$FID = predicted_expr$IID = NULL
  
  predicted_expr = predicted_expr[which(rownames(predicted_expr) %in% samples),]
  
  # Get performance predictions
  perf_epixcan = get(load(perf_file))
  
  # Subset predictions
  predicted_expr = predicted_expr[, which(colnames(predicted_expr) %in% perf_epixcan$gene)]
  
  return(list(
    predictions = predicted_expr,
    performance = perf_epixcan
  ))
}

#' Process INGENE Data
#' @description Process and prepare INGENE prediction data
#' @param dataset Name of the dataset
#' @param input.dir Input directory path
#' @param brain.region Brain region name
#' @return List containing processed INGENE data
process_ingene_data = function(dataset, input.dir, brain.region, samples) {
  # Construct file paths
  pred_file = file.path(input.dir, "INGENE/INGENE_averaged",
                       get_method_file_path("ingene", "predictions",
                                          brain_region = brain.region))
  perf_file = file.path("./predictions/LIBD_training/GTeX.v9/INGENE/INGENE_averaged",
                       get_method_file_path("ingene", "performance",
                                          dataset = "GTEx",
                                          brain_region = brain.region))
  
  # Read and process data
  predicted_expr = fread(pred_file) %>% as.data.frame()
  predicted_expr[[1]] = NULL
  rownames(predicted_expr) = predicted_expr$FID
  predicted_expr$FID = predicted_expr$IID = NULL
  
  predicted_expr = predicted_expr[which(rownames(predicted_expr) %in% samples),]
  
  # Get performance predictions
  perf_ingene = get(load(perf_file))
  
  # Subset predictions
  predicted_expr = predicted_expr[, which(colnames(predicted_expr) %in% perf_ingene$gene[perf_ingene$network == "Averaged_Ensemble"])]
  
  return(list(
    predictions = predicted_expr,
    performance = perf_ingene
  ))
}

#' Process CIS Data
#' @description Process and prepare CIS prediction data
#' @param dataset Name of the dataset
#' @param input.dir Input directory path
#' @param brain.region Brain region name
#' @return List containing processed CIS data
process_cis_data = function(dataset, input.dir, brain.region,samples) {
  # Construct file paths
  pred_file = file.path(input.dir, "CIS",
                       get_method_file_path("cis", "predictions",
                                          brain_region = brain.region,
                                          dataset = dataset))
  perf_file = file.path( "./predictions/LIBD_training/GTeX.v9/CIS",
                       get_method_file_path("cis", "performance",
                                          dataset = "GTEx",
                                          brain_region = brain.region))
  
  # Read and process data
  predicted_expr = fread(pred_file) %>% as.data.frame()
  rownames(predicted_expr) = predicted_expr$FID
  predicted_expr[[1]] = NULL
  predicted_expr$FID = predicted_expr$IID = NULL
  
  predicted_expr = predicted_expr[which(rownames(predicted_expr) %in% samples),]
  
  # Get performance predictions
  perf_cis = get(load(perf_file))
  
  # Subset predictions
  predicted_expr = predicted_expr[, which(colnames(predicted_expr) %in% perf_cis$gene)]
  
  return(list(
    predictions = predicted_expr,
    performance = perf_cis
  ))
}

###########################
# Gene Analysis Utils     #
###########################

#' Find Unique Genes Across Methods
#' @description Identify genes unique to each prediction method
#' @param module_data MODULE prediction data
#' @param epixcan_data EpiXcan prediction data
#' @param ingene_data INGENE prediction data
#' @param cis_data CIS prediction data
#' @return List of unique genes for each method
find_unique_genes = function(module_data, epixcan_data, ingene_data, cis_data) {
  # Extract gene lists
  module_genes = colnames(module_data$predictions)
  epixcan_genes = colnames(epixcan_data$predictions)
  ingene_genes = colnames(ingene_data$predictions)
  cis_genes = colnames(cis_data$predictions)
  
  # Find unique genes for each method
  unique_genes = list(
    module = setdiff(module_genes, union(union(epixcan_genes, ingene_genes), cis_genes)),
    epixcan = setdiff(epixcan_genes, union(union(module_genes, ingene_genes), cis_genes)),
    ingene = setdiff(ingene_genes, union(union(module_genes, epixcan_genes), cis_genes)),
    cis = setdiff(cis_genes, union(union(module_genes, epixcan_genes), ingene_genes))
  )
  
  return(unique_genes)
}

#' Find Common Genes Across Methods
#' @description Identify genes common to multiple prediction methods
#' @param module_data MODULE prediction data
#' @param epixcan_data EpiXcan prediction data
#' @param ingene_data INGENE prediction data
#' @param cis_data CIS prediction data
#' @return List of common genes for different method combinations
find_common_genes = function(module_data, epixcan_data, ingene_data, cis_data) {
  # Extract gene lists
  module_genes = colnames(module_data$predictions)
  epixcan_genes = colnames(epixcan_data$predictions)
  ingene_genes = colnames(ingene_data$predictions)
  cis_genes = colnames(cis_data$predictions)
  
  # Find common genes for different combinations
  common_genes = list(
    EMI = Reduce(intersect, list(epixcan_genes, module_genes, ingene_genes)),
    MIC = Reduce(intersect, list(module_genes, ingene_genes, cis_genes)),
    EM = setdiff(intersect(epixcan_genes, module_genes), ingene_genes),
    EI = setdiff(intersect(epixcan_genes, ingene_genes), module_genes),
    MI = setdiff(intersect(module_genes, ingene_genes), epixcan_genes),
    MC = setdiff(intersect(module_genes, cis_genes), ingene_genes),
    IC = setdiff(intersect(ingene_genes, cis_genes), module_genes)
  )
  
  return(common_genes)
}

#' Prepare Final Predictions
#' @description Combine predictions from different methods based on unique and common genes
#' @param unique_genes List of unique genes
#' @param common_genes List of common genes
#' @param module_data MODULE prediction data
#' @param epixcan_data EpiXcan prediction data
#' @param ingene_data INGENE prediction data
#' @param cis_data CIS prediction data
#' @return List of combined predictions
prepare_predictions = function(unique_genes, common_genes, module_data, epixcan_data,
                             ingene_data, cis_data) {
  # Prepare unique predictions
  unique_predictions = list(
    MODULE = module_data$predictions[, unique_genes$module, drop=FALSE],
    EpiXcan = epixcan_data$predictions[, unique_genes$epixcan, drop=FALSE],
    INGENE = ingene_data$predictions[, unique_genes$ingene, drop=FALSE],
    CIS = cis_data$predictions[, unique_genes$cis, drop=FALSE]
  )
  
  # Prepare common predictions
  common_predictions = list(
    EMI = list(
      MODULE = module_data$predictions[, common_genes$EMI, drop=FALSE],
      EpiXcan = epixcan_data$predictions[, common_genes$EMI, drop=FALSE],
      INGENE = ingene_data$predictions[, common_genes$EMI, drop=FALSE]
    ),
    MIC = list(
      MODULE = module_data$predictions[, common_genes$MIC, drop=FALSE],
      INGENE = ingene_data$predictions[, common_genes$MIC, drop=FALSE],
      CIS = cis_data$predictions[, common_genes$MIC, drop=FALSE]
    ),
    EM = list(
      MODULE = module_data$predictions[, common_genes$EM, drop=FALSE],
      EpiXcan = epixcan_data$predictions[, common_genes$EM, drop=FALSE]
    ),
    EI = list(
      EpiXcan = epixcan_data$predictions[, common_genes$EI, drop=FALSE],
      INGENE = ingene_data$predictions[, common_genes$EI, drop=FALSE]
    ),
    MI = list(
      MODULE = module_data$predictions[, common_genes$MI, drop=FALSE],
      INGENE = ingene_data$predictions[, common_genes$MI, drop=FALSE]
    ),
    MC = list(
      MODULE = module_data$predictions[, common_genes$MC, drop=FALSE],
      CIS = cis_data$predictions[, common_genes$MC, drop=FALSE]
    ),
    IC = list(
      INGENE = ingene_data$predictions[, common_genes$IC, drop=FALSE],
      CIS = cis_data$predictions[, common_genes$IC, drop=FALSE]
    )
  )
  
  return(list(
    unique_predictions = unique_predictions,
    common_predictions = common_predictions
  ))
}

###########################
# Main Analysis Function  #
###########################

#' Enhanced Data Preparation for Linear Model Analysis
#' @description Improved version of getDataforANOVA with better error handling and modularity
#' @param dataset Name of the dataset
#' @param input.dir Input directory path
#' @param true_expr True expression data
#' @param threshold Threshold for filtering
#' @param brain.region Brain region name
#' @return List containing processed data for enhanced analysis
getDataforLM = function(dataset, input.dir, true_expr, threshold, brain.region) {
  # Input validation
  if (is.null(dataset) || is.null(input.dir) || is.null(brain.region)) {
    stop("Required parameters cannot be NULL")
  }
  
  # Process data from each method
  print('Processing MODULE data...')
  module_data = process_module_data(dataset, input.dir, brain.region, samples=rownames(true_expr))
  samples = module_data$samples
  module_data$samples=NULL
  
  if(brain.region != "amygdala"){
    print('Processing EpiXcan data...')
    epixcan_data = process_epixcan_data(dataset, input.dir, brain.region, samples=samples)
  }else{epixcan_data=list()}
 
  print('Processing INGENE data...')
  ingene_data = process_ingene_data(dataset, input.dir, brain.region, samples=samples)
  
  print('Processing CIS data...')
  cis_data = process_cis_data(dataset, input.dir, brain.region, samples=samples)
  
  # Find unique and common genes
  unique_genes = find_unique_genes(module_data, epixcan_data, ingene_data, cis_data)
  common_genes = find_common_genes(module_data, epixcan_data, ingene_data, cis_data)
  
  # Prepare predictions
  predictions = prepare_predictions(unique_genes, common_genes, module_data, epixcan_data,
                                 ingene_data, cis_data)
  
  return(predictions)
}

###########################
# Model Fitting Functions #
###########################

#' Fit Linear Models for Gene Expression Prediction
#' @param df.lm Data frame containing predictors and true expression
#' @param combo String indicating which combination of predictors to use
#' @param gene_cols Named list of column indices for each predictor type
#' @return List of fitted linear models
fit_linear_models = function(df.lm, combo, gene_cols) {
  # Get the true expression column name (should be the last column)
  true_expr_col <- colnames(df.lm)[ncol(df.lm)]
  
  with(gene_cols, {
    switch(combo,
           "EMI" = list(
             model.EMI = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, MODULE.gene, INGENE.gene, true_expr_col)]),
             model.EM = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, MODULE.gene, true_expr_col)]),
             model.MI = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(INGENE.gene, MODULE.gene, true_expr_col)]),
             model.EI = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(INGENE.gene, EpiXcan.gene, true_expr_col)])
           ),
           "EM" = list(
             model.EM = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, MODULE.gene, true_expr_col)])
           ),
           "EI" = list(
             model.EI = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, INGENE.gene, true_expr_col)])
           ),
           "MI" = list(
             model.MI = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(MODULE.gene, INGENE.gene, true_expr_col)])
           ),
           "EMIC" = list(
             model.EMIC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, MODULE.gene, INGENE.gene, CIS.gene, true_expr_col)]),
             model.MIC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(MODULE.gene, INGENE.gene, CIS.gene, true_expr_col)]),
             model.EIC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, INGENE.gene, CIS.gene, true_expr_col)]),
             model.EMI = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, MODULE.gene, INGENE.gene, true_expr_col)]),
             model.EMC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, MODULE.gene, CIS.gene, true_expr_col)]),
             model.EM = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, MODULE.gene, true_expr_col)]),
             model.MI = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(INGENE.gene, MODULE.gene, true_expr_col)]),
             model.EI = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(INGENE.gene, EpiXcan.gene, true_expr_col)]),
             model.EC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(CIS.gene, EpiXcan.gene, true_expr_col)]),
             model.MC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(CIS.gene, MODULE.gene, true_expr_col)]),
             model.IC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(CIS.gene, INGENE.gene, true_expr_col)])
           ),
           "MIC" = list(
             model.MIC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(MODULE.gene, INGENE.gene, CIS.gene, true_expr_col)]),
             model.MI = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(INGENE.gene, MODULE.gene, true_expr_col)]),
             model.MC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(CIS.gene, MODULE.gene, true_expr_col)]),
             model.IC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(CIS.gene, INGENE.gene, true_expr_col)])
           ),
           "MEC" = list(
             model.EMC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(MODULE.gene, EpiXcan.gene, CIS.gene, true_expr_col)]),
             model.EM = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, MODULE.gene, true_expr_col)]),
             model.MC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(CIS.gene, MODULE.gene, true_expr_col)]),
             model.EC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, CIS.gene, true_expr_col)])
           ),
           "IEC" = list(
             model.EIC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(INGENE.gene, EpiXcan.gene, CIS.gene, true_expr_col)]),
             model.EI = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, INGENE.gene, true_expr_col)]),
             model.IC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(CIS.gene, INGENE.gene, true_expr_col)]),
             model.EC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, CIS.gene, true_expr_col)])
           ),
           "MC" = list(
             model.MC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(MODULE.gene, CIS.gene, true_expr_col)])
           ),
           "EC" = list(
             model.EC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(EpiXcan.gene, CIS.gene, true_expr_col)])
           ),
           "IC" = list(
             model.IC = lm(as.formula(paste(true_expr_col, "~ .")), data = df.lm[, c(INGENE.gene, CIS.gene, true_expr_col)])
           )
    )
  })
}

