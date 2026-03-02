
# Description: Analyzes and compares cis-derived predicted vs observed gene expression in testing data


# Load required libraries
required_packages = c("purrr", "limma", "dplyr", "tidyr", "tidyverse","data.table")
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Set working directory and base directory
# setwd("path/to/project/root")

# Load helper functions
source("./STEP8.Validate_cis_predictions_testing_data/validate_predictions_helper_functions.R")

#' Load and Process Testing Data
#' @param data_path Path to the testing RNA seq data
#' @param tissues Vector of tissue names to analyze
#' @return List of expression data for each tissue

testing_prediction_folder = ""
dataset_testing         = ""
dataset_training        = ""

load_testing_data = function(data_path = './rnaseq/your_testing_rnaseq.RData', 
                            tissues    = c("dlpfc")) {
  testing_data = get(load(data_path))
  map(tissues %>% set_names(.,.), ~ {
    tissue = .x
   testing_data[[tissue]]$assays$expression
    
  })
}

#' Analyze Gene Expression Predictions
#' @param model_type Type of model ("EpiXcan" or "CIS")
#' @param regions Vector of brain regions to analyze
#' @param true_expr True expression data
#' @param config List of configuration parameters
#' @return List of performance metrics for each region
analyze_predictions = function(model_type, dataset_training, testing_prediction_folder, name_target, 
                              regions,
                              true_expr,
                              config = list()) {
  
  # Set default configuration
  config = list(
      base_path = "./predictions/%s/",
      file_pattern = "%s_%s_%s_predicted.txt.gz",
      output_pattern = "%s_%s_performance_%s.RData"
     )

  
  # Function to process a single region
  process_region = function(region, model = NULL, dataset_training, testing_prediction_folder, name_target) {
    message(paste("Processing",dataset_training, model, "predictions in", name_target,"for region:", region))
    
    # Construct file paths based on model type
    pred_expr_file_path = file.path(sprintf(config$base_path, dataset_training, testing_prediction_folder), model)
  
    # Construct file names based on model type
    pred_expr_file_name = sprintf(config$file_pattern, model, region,name_target)
    
    ## output file name
    output_file = file.path(pred_expr_file_path, sprintf(config$output_pattern, model, region, name_target))
    if(file.exists(output_file)){return(NULL)}
    
    # Check if prediction file exists
    pred_file_path = file.path(pred_expr_file_path, pred_expr_file_name)
    if (!file.exists(pred_file_path)) {
      warning(paste("Prediction file not found:", pred_file_path))
      return(NULL)
    }
    
    # Read and analyze predictions
    pred_expr_df = fread(file.path(pred_expr_file_path,pred_expr_file_name))
    pred_expr_df = as.data.frame(pred_expr_df)
    
    if(region %in% c("sACC","dACC")){real.expr = true_expr[["acc"]]}else{real.expr=true_expr[[region]]}
    
    testing_data_pred_perf = run_model_analysis(
      real.expr            = real.expr,
      cov_df               = NULL,
      predicted.expr       = pred_expr_df,
      threshold            = "yes"
    )
    
    
    # Save results
    save(testing_data_pred_perf, file = output_file)
    
    return(testing_data_pred_perf)
  }
  
  # Process each region
  walk(regions %>% set_names(.,.), function(region) {
      process_region(region, model_type, dataset_training, testing_prediction_folder,name_target)
    })
  
}

# Main execution
main = function() {
  # Load testing data
  testing_data_true_expr = load_testing_data()
  
  if(dataset_testing == "CMC"){
    testing_data_true_expr = map(names(testing_data_true_expr) %>% set_names(names(testing_data_true_expr)), function(region){
      rownames(testing_data_true_expr[[region]]) = paste0(rownames(testing_data_true_expr[[region]]),"_", rownames(testing_data_true_expr[[region]]))
      return(testing_data_true_expr[[region]])
    })
  }
  
  # Analyze EpiXcan predictions
  epixcan_results = analyze_predictions(
    model_type    = "EpiXcan", dataset_training = dataset_training, testing_prediction_folder = testing_prediction_folder, name_target = dataset_testing,
    regions       = c("sACC","dACC","dlpfc"),
    true_expr     = testing_data_true_expr
  )

  # Analyze CIS  predictions
  cis_results  = analyze_predictions(
    model_type = "CIS", dataset_training = dataset_training, testing_prediction_folder = testing_prediction_folder, name_target = dataset_testing,
    regions    =  c("sACC","dACC","dlpfc"),
    true_expr  = testing_data_true_expr
  )
  return(NULL)
}

# Run the analysis
results = main()

