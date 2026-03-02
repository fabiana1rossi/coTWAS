#####################################################################
# INGENE Data Preparation Helper Functions
#####################################################################

#' Prepare data for a single brain region
#' @description Processes and combines different data sources for one brain region:
#' 1. Evaluates performance of EpiXcan and CIS models on testing_data data
#' 2. Selects best performing model for each gene
#' 3. Creates combined training dataset using selected models
#' 
#' @param testing.real.data Real testing_data expression data
#' @param epixcan_testing.imputed.data EpiXcan predictions for testing_data
#' @param cis_testing.imputed.data CIS predictions for testing_data
#' @param epixcan_training_data.imputed.data EpiXcan predictions for training_data
#' @param cis_training_data.imputed.data CIS predictions for training_data
#' @param training.real training_data real expression data
#' @return List containing processed training data
prepare_region_data = function(testing.real.data, epixcan_testing.imputed.data,
                             cis_testing.imputed.data, epixcan_training.imputed.data,
                             cis_training.imputed.data, training.real) {
    
    if(!is.null(epixcan_testing.imputed.data)){  ## not all regions have EpiXcan because annnotations were missing.
      # Evaluate EpiXcan performance
      epi_perf = evaluate_model_performance(
        real.expr      = testing.real.data,
        predicted.expr = epixcan_testing.imputed.data,
        model          = "EpiXcan"
      )
      
      # Evaluate CIS performance
      cis_perf = evaluate_model_performance(
        real.expr      = testing.real.data,
        predicted.expr = cis_testing.imputed.data,
        model          = "CIS"
      )
      
      # Select best performing model for each gene
      model_selection = select_best_model(epi_perf, cis_perf)
      
      ## Subset performance dataframes with corresponding winning genes
      cis_genes  = cis_perf  %>% filter(!genes %in% model_selection$genes[model_selection$win_model=="EpiXcan"]) %>% pull(genes)
      epixcan_genes  = epi_perf %>% filter(!genes %in% model_selection$genes[model_selection$win_model=="CIS"]) %>% pull(genes)
      
      
      # Create combined training dataset
      training_data = create_training_dataset(
        model_selection = model_selection,
        epixcan_data    = epixcan_training.imputed.data,
        cis_data        = cis_training.imputed.data,
        epixcan_genes   = epixcan_genes,
        cis_genes       = cis_genes,
        training_real   = training.real
      )
    }else{
      # Evaluate CIS performance
      cis_perf = evaluate_model_performance(
        real.expr      = testing.real.data,
        predicted.expr = cis_testing.imputed.data,
        model          = "CIS"
      )
      
      ## Subset performance dataframes with corresponding winning genes
      cis_genes  = cis_perf  %>% pull(genes)
      
      # Match samples between real and predicted CIS data
      common_samples   = intersect(rownames(training.real),
                                   cis_training.imputed.data$IID)
  
      training_real      = training.real[common_samples,]
      training_predicted = cis_training.imputed.data[which(cis_training.imputed.data$IID %in% common_samples),which(colnames(cis_training.imputed.data) %in% c("IID",cis_genes))]
      rownames(training_predicted) = training_predicted$IID
      training_predicted$IID = training_predicted$FID = NULL
      
      training_data = list(
        training.real = training_real,
        training.predicted = training_predicted,
        cis_genes          = cis_genes
      )
    }
    
    
    return(training_data)
}


step9InputDF <- function(real.expr, cov_df, output_dir=NULL, predicted.expr=NULL, predicted_df_file_path = NULL, 
                         samples=NULL,Age=NULL, dataset =NULL, threshold=NULL){
  
  
  ## Load predicted df. 
  if(!is.null(predicted_df_file_path)){
    predicted.expr =  read.table(predicted_df_file_path, header = T)
    if(!is.null(samples)){
      predicted.expr = predicted.expr[which(predicted.expr$FID %in% samples),]
    }
  }
  
  ids = predicted.expr$IID = gsub('\\.','-',predicted.expr$IID)
  predicted.expr[predicted.expr=='Inf'] = NA
  predicted.expr = as.data.frame(t(na.omit(t(predicted.expr))))
  predicted.expr = apply(predicted.expr,2,as.numeric)
  predicted.expr = as.data.frame(t(na.omit(t(predicted.expr))))
  predicted.expr$IID = predicted.expr$FID = ids
  print(dim(predicted.expr))
  
  # ## If expression has not been adjusted for covariates during pipeline...
  if(is.null(cov_df)){
    #print("already adjusted")
    if(sum(duplicated( limma::strsplit2(colnames(real.expr),'\\.')[,1])) > 0){
      expression = real.expr[,-which(duplicated( limma::strsplit2(colnames(real.expr),'\\.')[,1]))]
    }else{  expression = real.expr}
    colnames(expression) =  limma::strsplit2(colnames(expression),'\\.')[,1]
    genes = limma::strsplit2(colnames(predicted.expr),'\\.')[,1]
    expression = expression[,which(colnames(expression) %in% genes),drop=F]
    
    
  }else{
    exp = real.expr
    genes = limma::strsplit2(colnames(predicted.expr),'\\.')[,1]
    exp = exp[,which(colnames(exp) %in% genes)]
    # Adjust
    print('Adjusting...')
    expression = exp
    for (i in 1:length(colnames(exp))) {
      fit = lm(exp[,i] ~ cov_df)
      expression[,i] <- fit$residuals
    }
  }
  
  expression = as.data.frame(t(expression))
  print(dim(expression))
  
  colnames(expression) = gsub('\\.','-',colnames(expression))
  
  ######################
  # 1.Run LM analysis ##  
  ######################
  ## Make list with predicted df of the model/s used. 
  predicted_datasets <- list(predicted.expr)
  model_names <- list("Comparing Gene Expression")
  
  lm_results_list <- sapply(1:length(predicted_datasets), function(x){run_LM_analysis(expr.predicted_df = predicted_datasets[[x]],
                                                                                      expr.observed_df =  expression, 
                                                                                      model_name = model_names[[x]],
                                                                                      threshold = threshold)},simplify = F,USE.NAMES = T)
  # if(lm_results_list[[1]] == 'error'){return('error')}
  names(lm_results_list) = model_names
  
  ######################
  # 2.Store results   ##
  ######################
  
  ## Create a table to store the results. 
  result_table <- as.data.frame(t(sapply(lm_results_list, function(x){  
    mean = round(mean(x$r_squared_gene$r_square_adjusted),digits = 4)                                 ## r-sq adjusted
    mean2 = round(mean(x$r_squared_gene$r_squared),digits = 4)                                        ## r-sq 
    median = round(median(x$r_squared_gene$r_squared), digits = 4)
    max = round(max(x$r_squared_gene$r_squared),digits = 4)
    sd  = round(sd(x$r_squared_gene$r_squared),digits = 4)
    
    data.frame(mean_r2_adj = mean,  median_r2= median, mean_r2 = mean2, max_r2 = max, sd_r2 = sd)
  })))
  
  
  ## Create a table to store the results. 
  result_table_2 <- as.data.frame(t(sapply(lm_results_list, function(x){  
    adj_rsq = mean(x$r_squared_gene$r_square_adjusted)
    r_squared = mean(x$r_squared_gene$r_squared)
    rmse = mean(x$r_squared_gene$rmse)
    mae = mean(x$r_squared_gene$mae)
    
    
    data.frame(adj_rsq = adj_rsq,r_sq = r_squared,  rmse = rmse, mae = mae)
  })))
  
  # Plot max R-sq gene.
  plot_list <- sapply(1:length(lm_results_list), function(x){list(lm_results_list[[x]]$plot_best_gene)}, simplify = T, USE.NAMES = T)
  #do.call("grid.arrange", c(plot_list, ncol=2))
  # Plot R-sq distribution.
  plot_list_hist <- sapply(1:length(lm_results_list), function(x){list(lm_results_list[[x]]$hist)}, simplify = T, USE.NAMES = T)
  #do.call("grid.arrange", c(plot_list_hist, ncol=2))
  
  return(list(result_table, plot_list_hist, plot_list, as.data.frame(lm_results_list[[1]]$r_squared_gene), result_table_2))
}


#########################################################################
##  This function prepares datasets to be merged into a unique dataset ##
##  (genes as columns; participants as rows; label predicted/real).    ##
#########################################################################
#' @expr.predicted_df is a df with with the two first cols are called 
#' IID and FID respectively with sample IDs. The other cols are genes.
#' @expr.observed_df is the df obtained by the pre-processing step in the
#' "geneExprPreprocessR2" script. It has SampleIDs as columns and genes
#' as rows.
#' 

prepare_datasets = function(expr.predicted_df, expr.observed_df){
  expr.predicted_df = expr.predicted_df %>%
    dplyr::filter(., .$IID %in% colnames(expr.observed_df)) %>%
    dplyr::mutate(dataset = 'predicted')%>%
    dplyr::rename(.,participants = 'IID')%>%
    dplyr::select(-FID)
  
  if(sum(duplicated(limma::strsplit2(colnames(expr.predicted_df),'\\.')[,1]))!=0){
    expr.predicted_df = expr.predicted_df[,-which(duplicated(limma::strsplit2(colnames(expr.predicted_df),'\\.')[,1]))] 
  }
  colnames(expr.predicted_df) = limma::strsplit2(colnames(expr.predicted_df),"\\.")[,1]
  
  expr.observed_df = expr.observed_df[,which(colnames(expr.observed_df) %in% expr.predicted_df$participants)] %>%
    as.data.frame() %>%
    rownames_to_column(var = "genes") %>%
    pivot_longer(cols = -genes, names_to = "participants", values_to = "gene_expression")
  
  expr.observed_df = expr.observed_df %>%
    pivot_wider(names_from = genes, values_from = gene_expression) %>%
    dplyr::mutate(dataset = 'real') 
  
  # Dataset with ALL genes having expression and predicted expression for all participants. 
  merged_df = expr.predicted_df %>%
    bind_rows(expr.observed_df) 
  # Remove genes with all NAs in one or another dataset.
  genes_to_remove = c()
  for (i in 1:ncol(merged_df)){
    if(sum(is.na(merged_df[,i])) >= (sum(table(merged_df$dataset)/2))){
      genes_to_remove = append(genes_to_remove, i)
    }
  }
  if(!is.null(genes_to_remove)){
    merged_df = merged_df[,-genes_to_remove]
  }
  
  print(dim(merged_df))
  return(merged_df)
} 

###################################################################
## This function checks if there are any NA values.              ##
## If there are, subset the dataset with the non-missing values. ##
###################################################################
test_NA = function(merged_df){
  # How many genes w/o NA.
  test_NA_all_genes = merged_df %>%
    dplyr::select(-c(participants, dataset)) %>%
    map(~sum(!is.na(.))) %>%
    unlist() %>%
    as.data.frame() %>%
    rownames_to_column(var = 'gene_names') %>%
    dplyr::rename(.,non_missing = '.') 
  return(test_NA_all_genes)
}

#################################################################################################
##  Main function to: 1) define which genes are in both datasets and group them;               ##
##  2) for each nested group ( each gene has one col 'predicted'and one col 'real') perform    ##
##  linear analysis and return for each gene the adjusted R^2 value.                           ##
#################################################################################################
run_LM_analysis = function(expr.predicted_df, expr.observed_df, model_name, threshold=NULL){
  print(model_name)
  returnlist = list()
  merged_df = prepare_datasets(expr.predicted_df, expr.observed_df)
  test_NA_all_genes = test_NA(merged_df)
  
  
  ## Run the LM model.
  # Select relevant columns
  selected_cols = c("participants", "dataset", test_NA_all_genes$gene_names)
  summary_model_prediction = merged_df[, selected_cols]
  
  summary_model_prediction = merged_df %>%
    dplyr::select(participants, dataset, test_NA_all_genes$gene_names)%>%
    pivot_longer(cols = -c(participants, dataset), names_to = 'genes', values_to = 'expression') %>%
    pivot_wider(names_from = dataset, values_from = expression) %>%
    dplyr::group_nest(genes) %>%
    mutate(lm_test = map(data, ~(summary(lm(scale(real)~predicted, data = .x)))))
  
  
  r_squared_gene = summary_model_prediction %>%
    dplyr::select(genes)
  for(i in 1:nrow(summary_model_prediction)){
    r_squared_gene[i,2] = summary_model_prediction$lm_test[[i]]$adj.r.squared
    r_squared_gene[i,3] = summary_model_prediction$lm_test[[i]]$r.squared
    r_squared_gene[i,4] = summary_model_prediction$lm_test[[i]]$coefficients[2]
    r_squared_gene[i,5] = sqrt(mean(summary_model_prediction$lm_test[[i]]$residuals^2))
    r_squared_gene[i,6] = mean(abs(summary_model_prediction$lm_test[[i]]$residuals))
    r_squared_gene[i,7] = cor(summary_model_prediction$data[[i]]$real,summary_model_prediction$data[[i]]$predicted)
  }
  
  r_squared_gene = r_squared_gene %>%
    dplyr::rename(r_square_adjusted = '...2') %>%
    dplyr::rename(r_squared = '...3') %>%
    dplyr::rename(slope = "...4") %>%
    dplyr::rename(rmse = "...5")  %>%
    dplyr::rename(mae = "...6") %>%
    dplyr::rename(cor = "...7") 
  
  if(!is.null(threshold)){
    ## Take only positive correlation, positive slope and positive adj rsq
    rmse.quantile = quantile(r_squared_gene$rmse,0.75)
    
    r_squared_gene = r_squared_gene %>%
      dplyr::filter(., cor > 0,r_square_adjusted>0)
    print(paste0('# ', nrow(r_squared_gene), ' genes have positive correlation & adjR2!'))
    # r_squared_gene = r_squared_gene %>%
    #    dplyr::filter(.,  cor > 0)
    # print(paste0('# ', nrow(r_squared_gene), ' genes have positive correlation!'))
    if(nrow(r_squared_gene) == 0){return('error')}
    
  }else{print(paste0('# ', nrow(r_squared_gene), ' genes predicted!'))}
  
  r_squared_vec = c()
  summary_model_prediction = summary_model_prediction %>% filter(., genes %in% r_squared_gene$genes) 
  for(i in 1:nrow(summary_model_prediction)){
    r_squared_vec[i] = summary_model_prediction$lm_test[[i]]$r.squared
  }
  returnlist[["r_vec"]] =r_squared_vec
  
  # Plot R-sq distribution.
  r_squared_df = data.frame(r_squared_vec)
  print(round(max(r_squared_df),digits = 4))
  colnames(r_squared_df) = 'values'
  returnlist[["hist"]] = ggplot(r_squared_df, aes(x=values)) + 
    geom_histogram(binwidth=30) + 
    geom_histogram(color="white", fill="black") +
    xlim(0,max(r_squared_df)+0.1) + 
    ylim(0,8) +
    geom_vline(aes(xintercept=mean(values)),
               color="blue", linetype="dashed", size=0.7) +
    ggtitle(paste0(model_name, ': R-square values') )
  
  ## Look at the max R^2 and plot predicted vs real of the "best" predicted corresponding gene.
  best_gene = r_squared_gene[r_squared_gene$r_squared == max(r_squared_gene$r_squared),1]
  index_best_gene = which(summary_model_prediction[[1]] == unlist(best_gene))
  plot_db = as.data.frame(summary_model_prediction[[2]][[index_best_gene]]) %>%
    dplyr::mutate(residuals = real-predicted)
  
  returnlist[["plot_best_gene"]] = ggplot(plot_db, 
                                          aes(x = predicted,
                                              y = real)) + 
    geom_point() +
    geom_smooth(method = 'lm', color = 'black')+
    theme_classic(base_size = 10) +
    ggtitle(paste0(model_name,": Gene with highest r-squared"))
  
  returnlist[["r_squared_gene"]] = r_squared_gene
  
  return(returnlist)
}



#' Evaluate model performance
#' @param real.expr Real expression data
#' @param predicted.expr Predicted expression data
#' @param model Model name
#' @return Performance metrics data frame
evaluate_model_performance = function(real.expr, predicted.expr, model) {
    perf = step9InputDF(
        real.expr      = real.expr,
        predicted.expr = predicted.expr,
        dataset        = "testing_data",
        threshold      = 'yes', cov_df = NULL
    )[[4]]
    
    perf$Model = model
    return(perf)
}

#' Select best performing model for each gene
#' @param epi_perf EpiXcan performance metrics
#' @param cis_perf CIS performance metrics
#' @return Data frame with model selections
select_best_model = function(epi_perf, cis_perf) {
    model_selection = inner_join(epi_perf, cis_perf, by = "genes") %>%
        dplyr::rename(
            r_square_adjusted_epi = "r_square_adjusted.x",
            r_square_adjusted_cis = "r_square_adjusted.y"
        ) %>%
        mutate(win_model = case_when(
            r_square_adjusted_epi > r_square_adjusted_cis ~ "EpiXcan",
            r_square_adjusted_cis > r_square_adjusted_epi ~ "CIS",
            r_square_adjusted_epi == r_square_adjusted_cis ~ "Same"
        ))
    
    return(model_selection)
}

#' Create combined training dataset
#' @param model_selection Model selection results
#' @param epixcan_data EpiXcan predictions
#' @param cis_data CIS predictions
#' @param training_real Real training data
#' @return List containing processed training data
create_training_dataset = function(model_selection, epixcan_data, cis_data, epixcan_genes,
                                   cis_genes, training_real) {
    
    # Process EpiXcan data
    epixcan_subset = process_prediction_data(epixcan_data, epixcan_genes)
    
    # Process CIS data
    cis_subset    = process_prediction_data(cis_data, cis_genes)
    
    # Combine predictions
    training_predicted = combine_predictions(epixcan_subset, cis_subset)
    
    # Match samples between real and predicted data
    common_samples   = intersect(rownames(training_real),
                             rownames(training_predicted))

    training_real      = training_real[common_samples,]
    training_predicted = training_predicted[common_samples,]
    
    return(list(
        training.real = training_real,
        training.predicted = training_predicted,
        epixcan_genes      = epixcan_genes,
        cis_genes          = cis_genes
    ))
}

#' Process prediction data
#' @param data Prediction data
#' @param genes Genes to include
#' @return Processed prediction data
process_prediction_data = function(data, genes) {
    subset = data[, which(colnames(data) %in% c("FID", "IID", genes))]
    rownames(subset) = subset$FID
    return(subset)
}

#' Combine predictions from different models
#' @param epixcan_data EpiXcan predictions
#' @param cis_data CIS predictions
#' @return Combined prediction data
combine_predictions = function(epixcan_data, cis_data) {
    # Remove ID columns
    epixcan_data$FID = epixcan_data$IID = NULL
    cis_data$FID = cis_data$IID = NULL
    
    # Combine predictions
    combined = cbind(epixcan_data, cis_data)
    
    # Add back ID columns
    combined$FID = combined$IID = rownames(combined)
    
    return(combined)
} 