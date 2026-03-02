
#' Prepare and merge predicted and observed expression datasets
#' @param expr.predicted_df Dataframe containing predicted expression values
#' @param expr.observed_df Dataframe containing observed expression values
#' @return Merged dataframe with both predicted and observed values
prepare_datasets = function(expr.predicted_df, expr.observed_df) {
  # Filter and prepare predicted expression data
  expr.predicted_df = expr.predicted_df %>%
    dplyr::filter(IID %in% rownames(expr.observed_df)) %>% # Removed unnecessary . references
    dplyr::mutate(dataset = 'predicted') %>%
    dplyr::rename(participants = 'IID') %>% # Removed unnecessary . reference
    dplyr::select(-FID) %>%
    pivot_longer(cols = -c(participants, dataset), names_to = 'genes', values_to = 'gene_expression') # Simplified transformation
  
   # Reshape predicted data 
  expr.predicted_df = expr.predicted_df %>%
    pivot_wider(names_from = genes, values_from = gene_expression) %>%
    dplyr::mutate(dataset = 'predicted')
  
  ## Remove constant genes
  constant_elements = get_constant_columns(expr.predicted_df)
  expr.predicted_df = expr.predicted_df[,!(colnames(expr.predicted_df) %in% constant_elements)]
  
  # Prepare observed expression data
  expr.observed_df = expr.observed_df[which(rownames(expr.observed_df) %in% expr.predicted_df$participants),,drop=F] %>% t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "genes") %>%
    pivot_longer(cols = -genes, names_to = "participants", values_to = "gene_expression")
  
   # Reshape observed data to match predicted data format
  expr.observed_df = expr.observed_df %>% 
    pivot_wider(names_from = genes, values_from = gene_expression) %>%
    dplyr::mutate(dataset = 'real')
 
  ## Remove constant genes
  constant_elements = get_constant_columns(expr.observed_df)
  expr.observed_df = expr.observed_df[,!(colnames(expr.observed_df) %in% constant_elements)]
  
   
  # Merge datasets
  merged_df = expr.predicted_df %>%
    bind_rows(expr.observed_df) %>% as.data.frame()

  # Remove any rows containing NAs
 # merged_df = na.omit(merged_df)
  
 
  return(merged_df)
}


#' Test for NA values in gene expression data
#' @param merged_df Merged dataframe containing both predicted and observed expression values
#' @return List containing cleaned gene names and NA statistics
test_NA = function(merged_df) {
  # Get gene names (excluding participants and dataset columns)
  gene_names = colnames(merged_df)[!colnames(merged_df) %in% c("participants", "dataset")]
  
  # Calculate NA statistics for each gene
  na_stats = sapply(gene_names, function(gene) {
    gene_data = merged_df[[gene]]
    pred_data = gene_data[merged_df$dataset == "predicted"]
    obs_data = gene_data[merged_df$dataset == "real"]
    
    list(
      total_na = sum(is.na(gene_data)),
      pred_na = sum(is.na(pred_data)),
      obs_na = sum(is.na(obs_data)),
      prop_na = mean(is.na(gene_data))
    )
  })
  
  # Convert to data frame for easier handling
  na_stats_df = data.frame(
    gene = gene_names,
    total_na = unlist(na_stats["total_na",]),
    pred_na = unlist(na_stats["pred_na",]),
    obs_na = unlist(na_stats["obs_na",]),
    prop_na = unlist(na_stats["prop_na",])
  )
  
  # Print summary of NA statistics
  message("\nNA Statistics Summary:")
  message("Total genes checked: ", nrow(na_stats_df))
  message("Genes with any NA: ", sum(na_stats_df$total_na > 0))
  message("Genes with NA in predicted values: ", sum(na_stats_df$pred_na > 0))
  message("Genes with NA in observed values: ", sum(na_stats_df$obs_na > 0))
  
  # Return results
  return(list(
    gene_names = gene_names,
    na_stats = na_stats_df,
    genes_with_na = na_stats_df$gene[na_stats_df$total_na > 0],
    clean_genes = na_stats_df$gene[na_stats_df$total_na == 0]
  ))
} 



run_LM_analysis = function(expr.predicted_df, expr.observed_df, threshold=NULL) {
    
  # Initialize return list to store results
  results = list()
  
  # Prepare and clean datasets
  merged_df = prepare_datasets(expr.predicted_df, expr.observed_df)
  
  # Test for NA values across genes
  test_NA_all_genes = test_NA(merged_df)
  
  # Perform linear model analysis for each gene
  summary_model_prediction = run_gene_linear_models(merged_df, test_NA_all_genes$clean_genes)
  
  # Calculate performance metrics for each gene
  r_squared_gene = calculate_gene_metrics(summary_model_prediction)
  r_squared_gene = as.data.frame(r_squared_gene)

  # Apply threshold filtering if specified
  if (!is.null(threshold)) {
    r_squared_gene = filter_by_threshold(r_squared_gene)
    if (nrow(r_squared_gene) == 0) return(NULL)
  }
  
  return(r_squared_gene)
}

# Helper functions

get_constant_columns = function(df) {
  constant_elements = function(x) length(unique(x)) == 1
  constant_cols     = sapply(df[df$dataset == 'predicted',], constant_elements)
  names(constant_cols[constant_cols])[!grepl('dataset', names(constant_cols[constant_cols]))]
}

run_gene_linear_models = function(merged_df, gene_names) {
  merged_df %>%
    dplyr::select(participants, dataset, gene_names) %>%
    pivot_longer(cols = -c(participants, dataset), names_to = 'genes', values_to = 'expression') %>%
    pivot_wider(names_from = dataset, values_from = expression) %>%
    dplyr::group_nest(genes) %>%
    mutate(lm_test = map(data, ~(summary(lm(scale(real)~scale(predicted), data = .x)))))
}

calculate_gene_metrics = function(summary_model_prediction) {

 # Define NRMSE function
  NRMSE = function(residuals, actual) {
    rmse      = sqrt(mean(residuals^2, na.rm = TRUE))
    sd_actual = sd(actual, na.rm = TRUE)
    return(rmse / sd_actual)
  }
  
  # Initialize matrix with 9 columns
  r_squared_gene = matrix(NA, nrow = nrow(summary_model_prediction), ncol = 7)
  
  for(i in 1:nrow(summary_model_prediction)) {
    model_summary = summary_model_prediction$lm_test[[i]]
    real_vals     = summary_model_prediction$data[[i]]$real
    pred_vals      = summary_model_prediction$data[[i]]$predicted

    r_squared_gene[i, 1:7] = c(
      summary_model_prediction$genes[[i]],
      model_summary$adj.r.squared,
      model_summary$r.squared,
      ifelse(nrow(model_summary$coefficients) > 1, model_summary$coefficients[2, 1], NA),
      NRMSE(model_summary$residuals, real_vals),
      mean(abs(model_summary$residuals), na.rm = TRUE),
      cor(real_vals, pred_vals, use = "complete.obs")
      
    )
  }
  
  # Assign column names
  colnames(r_squared_gene) = c(
    "gene", "r_square_adjusted", "r_squared", "slope",
    "nrmse", "mae", "cor"
  )
  
  # Convert numeric columns
  r_squared_gene[, 2:ncol(r_squared_gene)] = apply(r_squared_gene[, 2:ncol(r_squared_gene)], 2, as.numeric)
  
  return(as.data.frame(r_squared_gene))
}

filter_by_threshold = function(r_squared_gene) {
  
  filtered = r_squared_gene %>% 
    dplyr::filter(cor > 0, r_square_adjusted > 0)
  print(paste0('# ', nrow(filtered), ' genes have positive correlation & adj rsquared'))
  filtered
}

get_rsquared_vector = function(summary_model_prediction, r_squared_gene) {
  summary_model_prediction %>%
    filter(genes %in% r_squared_gene$genes) %>%
    pull(lm_test) %>%
    sapply(function(x) x$r.squared)
}

plot_rsquared_distribution = function(r_squared_vec, model_name) {
  r_squared_df = data.frame(values = r_squared_vec)
  print(round(mean(r_squared_vec), digits = 4))
  
  ggplot(r_squared_df, aes(x=values)) + 
    geom_histogram(color="white", fill="black", binwidth=30) +
    xlim(0, max(r_squared_df) + 0.1) + 
    ylim(0, 8) +
    geom_vline(aes(xintercept=mean(values)),
               color="blue", linetype="dashed", size=0.7) +
    ggtitle(paste0(model_name, ': R-square values'))
}

plot_best_predicted_gene = function(summary_model_prediction, r_squared_gene, model_name) {
  best_gene = r_squared_gene[r_squared_gene$r_squared == max(r_squared_gene$r_squared), 1]
  index_best_gene = which(summary_model_prediction[[1]] == unlist(best_gene))
  
  plot_db = as.data.frame(summary_model_prediction[[2]][[index_best_gene]]) %>%
    dplyr::mutate(residuals = real - predicted)
  
  ggplot(plot_db, aes(x = predicted, y = real)) + 
    geom_point() +
    geom_smooth(method = 'lm', color = 'black') +
    theme_classic(base_size = 10) +
    ggtitle(paste0(model_name, ": Gene with highest r-squared"))
}

run_model_analysis = function(real.expr, cov_df,  predicted.expr=NULL, threshold=NULL) {
  
  # Clean and process predicted expression data
  predicted.expr = process_predicted_expression(predicted.expr)
  
  # Process real expression data based on covariates
  expression = process_real_expression(real.expr, cov_df, predicted.expr)
  
  # Run linear model analysis
  lm_results =  run_LM_analysis(expr.predicted_df = predicted.expr,
                                expr.observed_df  = expression,
                                threshold         = threshold)
  
  # Generate and format results
  return(lm_results)
}


# Helper functions for run_model_analysis
load_predicted_expression = function(predicted_df_file_path, predicted.expr, extra_df) {
  # Load from file if path provided
  if (!is.null(predicted_df_file_path)) {
    predicted.expr = read.table(predicted_df_file_path, header = T)
  }
  
  # Filter genes based on extra_df if provided
  if (!is.null(extra_df)) {
    predicted.expr = predicted.expr[, which(colnames(predicted.expr) %in% 
                                             c("FID", "IID", extra_df$gene)), 
                                   drop = FALSE]
  }
  
  return(predicted.expr)
}

process_predicted_expression = function(predicted.expr) {
  # Extract and preserve IDs
  ids = if (!is.null(predicted.expr$IID)) predicted.expr$IID else rownames(predicted.expr)
  
  # Clean infinite values and convert to numeric
  predicted.expr[predicted.expr == 'Inf'] = NA
  predicted.expr = as.data.frame(t(na.omit(t(predicted.expr))))
  predicted.expr = apply(predicted.expr, 2, as.numeric)
  predicted.expr = as.data.frame(t(na.omit(t(predicted.expr))))
  
  # Restore IDs
  predicted.expr$IID = predicted.expr$FID = ids
  
  print(paste("Processed predicted expression dimensions:", 
              paste(dim(predicted.expr), collapse=" x ")))
  
  return(predicted.expr)
}

process_real_expression = function(real.expr, cov_df, predicted.expr) {
  if (is.null(cov_df)) {
    # Handle pre-adjusted expression data
    expression = handle_preadjusted_expression(real.expr, predicted.expr)
  } else {
    # Adjust expression for covariates
    expression = adjust_for_covariates(real.expr, predicted.expr, cov_df)
  }
  
  print(paste("Processed real expression dimensions:", 
              paste(dim(expression), collapse=" x ")))
  
  return(expression)
}

handle_preadjusted_expression = function(real.expr, predicted.expr) {
  # Remove duplicates if present
  if (sum(duplicated(limma::strsplit2(colnames(real.expr), '\\.')[,1])) > 0) {
    expression = real.expr[, -which(duplicated(limma::strsplit2(colnames(real.expr), '\\.')[,1]))]
  } else {
    expression = real.expr
  }
  
  # Clean column names and filter genes
  colnames(expression) = limma::strsplit2(colnames(expression), '\\.')[,1]
  genes = limma::strsplit2(colnames(predicted.expr), '\\.')[,1]
  expression = expression[, which(colnames(expression) %in% genes), drop=FALSE]
  
  return(expression)
}

adjust_for_covariates = function(real.expr, predicted.expr, cov_df) {
  # Filter genes
  genes = limma::strsplit2(colnames(predicted.expr), '\\.')[,1]
  exp = real.expr[, which(colnames(real.expr) %in% genes)]
  
  # Adjust for covariates
  print('Adjusting for covariates...')
  expression = exp
  for (i in 1:length(colnames(exp))) {
    fit = lm(exp[,i] ~ cov_df)
    expression[,i] = fit$residuals
  }
  
  return(expression)
}



create_summary_table = function(lm_results_list) {
  as.data.frame(t(sapply(lm_results_list, function(x) {
    mean = round(mean(x$r_squared_gene$r_square_adjusted), digits = 4)
    mean2 = round(mean(x$r_squared_gene$r_squared), digits = 4)
    median = round(median(x$r_squared_gene$r_squared), digits = 4)
    max = round(max(x$r_squared_gene$r_squared), digits = 4)
    sd = round(sd(x$r_squared_gene$r_squared), digits = 4)
    
    data.frame(mean_r2_adj = mean, median_r2 = median, 
               mean_r2 = mean2, max_r2 = max, sd_r2 = sd)
  })))
}

create_detailed_table = function(lm_results_list) {
  as.data.frame(t(sapply(lm_results_list, function(x) {
    data.frame(
      adj_rsq = mean(x$r_squared_gene$r_square_adjusted),
      r_sq = mean(x$r_squared_gene$r_squared),
      rmse = mean(x$r_squared_gene$rmse),
      mae = mean(x$r_squared_gene$mae)
    )
  })))
}




