#####################################################################
# Helper Functions for Training Pipeline (Basic Elastic Net)
#####################################################################

# Helper function to get file paths
get_file_path = function(file_type, brain_region=NULL, chr=NULL, CONFIG) {
  pattern = CONFIG$file_patterns[[file_type]]
  
  if(file_type %in% c("geno")) {
    return(sprintf(file.path(CONFIG$genotype_by_chr_dir, brain_region, pattern), chr))
  } else if(file_type %in% c("eqtl")) {
    return(sprintf(file.path(CONFIG$eqtls_dir, brain_region, pattern), chr))
  } else if(file_type %in% c("weights")) {
    return(sprintf(file.path(CONFIG$output_dir,"weights", brain_region, pattern), 
                   brain_region, chr))
  } else if(file_type %in% c("summary")) {
    return(sprintf(file.path(CONFIG$output_dir, "summary", brain_region, pattern), 
                   brain_region, chr))
  }else if(file_type %in% c("cv_folds")) {
    return(sprintf(file.path(CONFIG$data_dir,pattern), brain_region))
  }
}



generate_fold_ids = function(n_samples, n_folds) {
  n = ceiling(n_samples / n_folds)
  fold_ids = rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}

adjust_for_covariates = function(expression, cov_df) {
  if(!is.null(cov_df)){
    print("Correcting expression")
    cov_df   = cov_df[,!grepl("BrNum|RNum|Region|Race",colnames(cov_df))]
    samples  = rownames(cov_df)
    cov_df$RIN = as.numeric(as.character(cov_df$RIN))
    cov_df$mitoRate = as.numeric(as.character(cov_df$mitoRate))
    cov_df$Age = as.numeric(as.character(cov_df$Age))
    cov_df$rRNA_rate = as.numeric(as.character(cov_df$rRNA_rate))
    cov_df$overallMapRate = as.numeric(as.character(cov_df$overallMapRate))
    cov_df$totalAssignedGene = as.numeric(as.character(cov_df$totalAssignedGene))
    cov_df$Protocol    = as.factor(cov_df$Protocol)
    cov_df$Dataset    = as.factor(cov_df$Dataset)
    
    cov_df[,11:ncol(cov_df)]   = apply(cov_df[,11:ncol(cov_df)],2,as.numeric)
    
    cov_df   = as.data.frame(lapply(cov_df, function(x) if (is.character(x)) as.factor(x) else x))
    rownames(cov_df) = samples
    
    stopifnot(identical(rownames(cov_df),rownames(expression)))
    
    combined_df = cbind(expression[,1], cov_df)
    colnames(combined_df)[1] = "expression_vec"
    
    expr_resid = as.data.frame(summary(lm(expression_vec ~ ., data=combined_df))$residuals)
    
  }else{
    expr_resid = expression[,1,drop=F]
  }
  
  expr_resid = scale(expr_resid, center = TRUE, scale = TRUE)
  expr_resid
}

calc_R2 = function(y, y_pred) {
  tss = sum(y**2)
  rss = sum((y - y_pred)**2)
  1 - rss/tss
}

calc_corr = function(y, y_pred) {
  sum(y*y_pred) / (sqrt(sum(y**2)) * sqrt(sum(y_pred**2)))
}

nested_cv_elastic_net_perf = function(x, y, n_samples, train_test_fold_ids, alpha, n_k_folds=4) {
  # Gets performance estimates for k-fold cross-validated elastic-net models.
  # Splits data into n_train_test_folds disjoint folds, roughly equal in size,
  # and for each fold, calculates a n_k_folds cross-validated elastic net model. Lambda parameter is
  # cross validated. Then get performance measures for how the model predicts on the hold-out
  # fold. Get the coefficient of determination, R^2, and a p-value, where the null hypothesis
  # is there is no correlation between prediction and observed.
  #
  # The mean and standard deviation of R^2 over all folds is then reported, and the p-values
  # are combined using Fisher's method.
  n_train_test_folds = length(unique(train_test_fold_ids))
  rmse_folds <- rep(0, n_train_test_folds)
  R2_folds <- rep(0, n_train_test_folds)
  corr_folds <- rep(0, n_train_test_folds)
  zscore_folds <- rep(0, n_train_test_folds)
  pval_folds <- rep(0, n_train_test_folds)
  # Outer-loop split into training and test set.
  y=as.vector(y)
  set.seed(070892)
  for (test_fold in 1:n_train_test_folds) {
    train_idxs <- which(train_test_fold_ids != test_fold)
    test_idxs <- which(train_test_fold_ids == test_fold)
    x_train <- x[train_idxs, ]
    y_train <- y[train_idxs]
    x_test <- x[test_idxs, ]
    y_test <- y[test_idxs]
    # Inner-loop - split up training set for cross-validation to choose lambda.
    # this is not related to the adjusted R2 between two models, preset cv ids are not passed  
    cv_fold_ids <- generate_fold_ids(length(y_train), n_k_folds)
    y_pred <- tryCatch({
      # Fit model with training data.
      # No penalty factors for basic elastic net
      fit <- cv.glmnet(x_train, y_train, nfolds = n_k_folds, alpha = alpha, type.measure='mse', foldid = cv_fold_ids)
      # Predict test data using model that had minimal mean-squared error in cross validation.
      predict(fit, x_test, s = 'lambda.min')},
      # if the elastic-net model did not converge, predict the mean of the y_train (same as all non-intercept coef=0)
      error = function(cond) rep(mean(y_train), length(y_test)))
    R2_folds[test_fold] <- calc_R2(y_test, y_pred)
    ## get rmse
    res <- summary(lm(y_test~y_pred))
    rmse_folds[test_fold] <-sqrt(1-res$r.squared)*(res$sigma)
    # Get p-value for correlation test between predicted y and actual y.
    # If there was no model, y_pred will have var=0, so cor.test will yield NA.
    # In that case, give a random number from uniform distribution, which is what would
    # usually happen under the null.
    corr_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor(y_pred, y_test), 0)
    zscore_folds[test_fold] <- atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3) # Fisher transformation
    pval_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
  }
  rmse_avg <- mean(rmse_folds)
  R2_avg <- mean(R2_folds)
  R2_sd <- sd(R2_folds)
  rho_avg <- mean(corr_folds)
  rho_se <- sd(corr_folds)
  rho_avg_squared <- rho_avg**2
  # Stouffer's method for combining z scores.
  zscore_est <- sum(zscore_folds) / sqrt(n_train_test_folds)
  zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
  # Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
  list(rmse_avg=rmse_avg,R2_avg=R2_avg, R2_sd=R2_sd, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval)
}

#' Filter genotype data to remove NAs
#' @param geno    genotype 
#' @return Filtered genotype matrix
filter_genotypes = function(geno) {
  geno = geno[, apply(geno, 2, function(x) all(!is.na(x)))]
  
  return(geno)
}

#' Get gene coordinates and annotation
#' @param gene_annot Gene annotation data frame
#' @param gene Gene identifier
#' @return List with gene coordinates and info
get_gene_info = function(gene_annot, gene) {
  gene_data = na.omit(gene_annot[gene_annot$gencodeID == gene,])
  if(nrow(gene_data) == 0) {
    warning(sprintf("No annotation found for gene %s", gene))
    return(list(
      coords = c(NA, NA),
      name = NA,
      type = NA
    ))
  }
  list(
    coords = c(gene_data$start, gene_data$end),
    name = gene_data$gene_name,
    type = gene_data$gene_type
  )
}

#' Get SNPs within cis window of a gene
#' @param genotypes Genotype matrix
#' @param snp_annot SNP annotation data frame
#' @param gene_coords Gene coordinates
#' @param cis_window Size of cis window
#' @return Matrix of cis SNPs
filter_cis_snps = function(genotypes, snp_annot, gene_coords, cis_window) {
  cis_snps = snp_annot %>%
    filter(
      Position >= (gene_coords[1] - cis_window),
      Position <= (gene_coords[2] + cis_window)
    ) %>%
    pull(SNP)
  
  if(length(cis_snps) == 0) return(NULL)
  
  return(genotypes[, which(colnames(genotypes) %in% cis_snps), drop=FALSE])
}

#' Process a single chromosome
#' @param brain_region Brain region name
#' @param chr Chromosome number
#' @param data_region List of input data
#' @param CONFIG Configuration parameters
process_chromosome = function(brain_region, chr, data_region, CONFIG) {
  message(sprintf("Processing chromosome %d for %s", chr, brain_region))
  
  # Skip if results exist
  if(file.exists(get_file_path("summary", brain_region, chr=chr, CONFIG = CONFIG))) {
    message("Results already exist, skipping...")
    return(NULL)
  }
  
  # Load chromosome-specific data
  genotypes_snp_annot = readRDS(get_file_path("geno", chr=chr, brain_region = brain_region, CONFIG = CONFIG))
  if(is.null(genotypes_snp_annot) || !all(c("genotype", "snp_info") %in% names(genotypes_snp_annot))) {
    stop(sprintf("Invalid genotype data format for chromosome %d", chr))
  }
  
  genotypes = genotypes_snp_annot$genotype
  snp_annot = genotypes_snp_annot$snp_info
  
  # Validate data structure
  if(is.null(genotypes) || is.null(snp_annot)) {
    stop(sprintf("Missing genotype or SNP annotation data for chromosome %d", chr))
  }
  
  if(nrow(snp_annot) != ncol(genotypes)) {
    stop(sprintf("Mismatch between number of SNPs in genotype (%d) and annotation (%d) for chromosome %d", 
                 ncol(genotypes), nrow(snp_annot), chr))
  }
  
  ## Remove SNP with NAs 
  filtered_genotypes = filter_genotypes(genotypes)
  if(ncol(filtered_genotypes) == 0) {
    warning(sprintf("No valid SNPs remaining after filtering for chromosome %d", chr))
    return(NULL)
  }
  
  ## Extract expression, covariates and genotype data
  if(is.null(data_region$assays) || is.null(data_region$assays$expression)) {
    stop("Missing expression data in data_region")
  }
  if(is.null(data_region$colData)) {
    stop("Missing covariate data in data_region")
  }
  if(is.null(data_region$rowData)) {
    stop("Missing gene annotation data in data_region")
  }
  
  expression = as.data.frame(data_region$assays$expression)
  covariates = as.data.frame(data_region$colData)
  gene_annot = as.data.frame(data_region$rowData)
  
  # Validate numeric columns
  gene_annot$chr   = as.numeric(gsub("chr","",gene_annot$chr))
  gene_annot       = na.omit(gene_annot) 
  gene_annot$start = as.numeric(gene_annot$start)
  gene_annot$end   = as.numeric(gene_annot$end)
  
  if(any(is.na(gene_annot$chr)) || any(is.na(gene_annot$start)) || any(is.na(gene_annot$end))) {
    warning("Some gene coordinates are NA - these genes will be skipped")
  }
  
  ## Match expression and covariates
  common_samples = intersect(rownames(filtered_genotypes), rownames(expression))
  if(length(common_samples) == 0) {
    stop("No common samples between genotype and expression data")
  }
  
  expression         = expression[common_samples, , drop=FALSE]
  filtered_genotypes = filtered_genotypes[common_samples, , drop=FALSE]
  covariates         = covariates[common_samples, , drop=FALSE]
  
  # Ensure consistent sample ordering
  covariates = covariates[match(common_samples, rownames(covariates)), , drop=FALSE]
  expression = expression[match(common_samples, rownames(expression)), , drop=FALSE]
  
  stopifnot(identical(rownames(expression), rownames(filtered_genotypes)))  
  stopifnot(identical(rownames(covariates), rownames(filtered_genotypes)))  
  
  # Load cis-eQTLs if available
  cis_eqtl_file       = get_file_path("eqtl", brain_region, chr, CONFIG)
  significant_ciseQTL = if(file.exists(cis_eqtl_file)) {
    readRDS(cis_eqtl_file)
  } else {
    NULL
  }

  ## Select cisEQTLs based on pvalue
  significant_ciseQTL2 = significant_ciseQTL %>% mutate(FDR = p.adjust(meta_p, method="fdr")) %>% filter(FDR < 0.05)

  # Process each gene
  genes = gene_annot$gencodeID[gene_annot$chr == chr]

  
  if(length(genes) == 0){
    warning(sprintf("No genes found for chromosome %d", chr))
    return(NULL)
  }
  
  message(sprintf("Processing %d genes for chromosome %d", length(genes), chr))
  
  ## Train Elastic Net 
  results = mclapply(genes, function(g) {
    message(sprintf("Processing gene: %s", g))
    tryCatch({
      result = train_gene_model(
        gene                = g,
        gene_annot          = gene_annot,
        expression          = expression,
        genotypes           = filtered_genotypes,
        covariate           = covariates,
        snp_annot           = snp_annot,
        ciseqtl             = significant_ciseQTL,
        brain_region        = brain_region,
        CONFIG              = CONFIG
      )
      if(is.null(result)) {
        message(sprintf("Gene %s returned NULL result", g))
      } else {
        message(sprintf("Gene %s processed successfully", g))
      }
      result
    }, error = function(e) {
      warning(sprintf("Error processing gene %s: %s", g, e$message))
      NULL
    })
  }, mc.cores = CONFIG$n_cores, mc.preschedule = F)
  
  # Ensure results is a list
  if (!is.list(results)) {
    warning("Results from mclapply is not a list, converting...")
    results = as.list(results)
  }

  
  # Create summary and weights files
  valid_results = lapply(results, function(gene_result) {
    # Check if gene_result is a list and has the expected structure
    if (!is.null(gene_result) && is.list(gene_result) && !is.null(gene_result$gene_info) && 
        !is.null(gene_result$model) && !is.null(gene_result$model$cv_results)) {
      data.frame(
        gene_id                = gene_result$gene,
        gene_name              = ifelse(is.na(gene_result$gene_info$name), "NA", gene_result$gene_info$name),
        gene_type              = ifelse(is.na(gene_result$gene_info$type), "NA", gene_result$gene_info$type),
        chromosome             = chr,
        alpha                  = CONFIG$alpha,
        n_snps_in_window       = gene_result$model$cv_results$n_snps_in_window,
        n_snps_in_model        = gene_result$model$cv_results$n_snps_in_model,
        best_lambda            = gene_result$model$cv_results$best_lambda,
        inner_R2_avg           = gene_result$model$cv_results$inner_R2_avg,
        inner_R2_sd            = gene_result$model$cv_results$inner_R2_sd,
        inner_pval_est         = gene_result$model$cv_results$inner_pval_est,
        inner_rho_avg          = gene_result$model$cv_results$inner_rho_avg,
        inner_rho_se           = gene_result$model$cv_results$inner_rho_se,
        inner_rho_zscore       = gene_result$model$cv_results$inner_rho_zscore,
        inner_rho_avg_squared  = gene_result$model$cv_results$inner_rho_avg_squared,
        inner_zscore_pval      = gene_result$model$cv_results$inner_zscore_pval,
        cv_R2_gene_avg         = gene_result$model$cv_results$cv_R2_gene_avg,
        cv_R2_gene_sd          = gene_result$model$cv_results$cv_R2_gene_sd,
        cv_rho_gene_avg        = gene_result$model$cv_results$cv_rho_gene_avg,
        training_R2            = gene_result$model$cv_results$training_R2,
        cv_rho_gene_se         = gene_result$model$cv_results$cv_rho_gene_se,
        cv_rho_gene_avg_squared = gene_result$model$cv_results$cv_rho_gene_avg_squared,
        cv_zscore_est_gene      = gene_result$model$cv_results$cv_zscore_est_gene,
        cv_zscore_pval_gene     = gene_result$model$cv_results$cv_zscore_pval_gene,
        cv_pval_est_gene        = gene_result$model$cv_results$cv_pval_est_gene,
        cor_gene_all_data_pred  = gene_result$model$cv_results$cor_gene_all_data_pred,
        adj_rsq_gene_all_data   = gene_result$model$cv_results$adj_rsq_gene_all_data,
        rmse_avg               = gene_result$model$cv_results$rmse_avg,
        pval_gene_all_data      = gene_result$model$cv_results$pval_gene_all_data,
        cv_adj_r2_gene_avg      = gene_result$model$cv_results$cv_adj_r2_gene_avg,
        cv_adj_r2_gene_sd       = gene_result$model$cv_results$cv_adj_r2_gene_sd,
        cv_pval_lm_gene_avg     = gene_result$model$cv_results$cv_pval_lm_gene_avg
      )
    } else {
      NULL
    }
  })
  
  # Remove NULL results and combine
  valid_results = valid_results[!sapply(valid_results, is.null)]
  
  if(length(valid_results) == 0) {
    warning(sprintf("No valid results found for chromosome %d", chr))
    model_summaries = data.frame()
  } else {
    model_summaries = do.call(rbind, valid_results)
  }
  
  # Save model summary
  summary_file = get_file_path("summary", brain_region, chr, CONFIG)
  
  write.table(
    model_summaries,
    file = summary_file,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  
  # Create weights file
  valid_weights = lapply(results, function(gene_result) {
    # Check if gene_result is a list and has the expected structure
    if (!is.null(gene_result) && is.list(gene_result) && !is.null(gene_result$model) && 
        !is.null(gene_result$model$weighted_snps_info) && nrow(gene_result$model$weighted_snps_info) > 0) {
      gene_result$model$weighted_snps_info
    } else {
      NULL
    }
  })
  
  # Remove NULL results and combine
  valid_weights = valid_weights[!sapply(valid_weights, is.null)]
  
  if(length(valid_weights) == 0) {
    warning(sprintf("No valid weights found for chromosome %d", chr))
    weights_data = data.frame()
  } else {
    weights_data = do.call(rbind, valid_weights)
  }
  
  # Save weights
  weights_file = get_file_path("weights", brain_region, chr, CONFIG)
 
  write.table(
    weights_data,
    file = weights_file,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
}

#' Train model for a single gene
#' @param gene Gene identifier
#' @param data List of input data
#' @param genotypes Genotype matrix
#' @param snp_annot SNP annotation
#' @param significant_ciseQTL Significant cis-eQTLs
train_gene_model = function(gene, gene_annot, expression, genotypes, covariate, snp_annot, ciseqtl, brain_region, CONFIG) {
  # Input validation
  if(is.null(gene) || is.null(gene_annot) || is.null(expression) || 
     is.null(genotypes) || is.null(snp_annot)) {
    stop("Missing required input data")
  }
  
  ## Get gene info
  gene_info = get_gene_info(gene_annot, gene)
  if(is.null(gene_info) || any(is.na(gene_info$coords))) {
    warning(sprintf("Invalid gene coordinates for gene %s", gene))
    return(NULL)
  }
  
  ################
  ## Get cis-SNPs
  start = gene_info$coords[1] - 1e6
  end   = gene_info$coords[2] + 1e6
  # Pull cis-SNP info
  cissnps = subset(snp_annot, snp_annot$Position >= start & snp_annot$Position <= end)
  
  # Pull cis-SNP genotypes
  cissnp1index = intersect(colnames(genotypes),cissnps$SNP)
  cis_geno <- genotypes[,cissnp1index, drop = FALSE]
  
  cm <- colMeans(cis_geno, na.rm = TRUE)
  minorsnps <- subset(colMeans(cis_geno), cm > 0 & cm < 2)
  minorsnps <- names(minorsnps)
  cis_geno <- cis_geno[,minorsnps, drop = FALSE]

  if(ncol(cis_geno) == 0) {
      warning(sprintf("No matching SNPs found in genotype data for gene %s", gene))
      return(NULL)
    }
  
  
  if(ncol(cis_geno) < 2) {
    warning(sprintf("Insufficient SNPs (< 2) for gene %s", gene))
    return(NULL)
  }

  ########### End cis-SNP selection
  
  # Train model
  model_results = tryCatch({
    train_elastic_net(
      gene_id         = gene,
      gene_annot      = gene_annot,
      expression      = expression,
      genotypes       = cis_geno,
      snp_annot       = snp_annot,
      covariates      = NULL,  ## Expression already adjusted for covs
      brain_region    = brain_region,
      CONFIG          = CONFIG
    )
  }, error = function(e) {
    warning(sprintf("Error training model for gene %s: %s", gene, e$message))
    NULL
  })
  
  if(is.null(model_results)) {
    return(NULL)
  }
  
  list(
    gene      = gene,
    gene_info = gene_info,
    model     = model_results
  )
}

#' Train elastic net model for a gene (basic version without priors)
#' @param gene_id Gene identifier
#' @param expression Expression data
#' @param genotypes Genotype data
#' @param snp_annot SNP annotation
#' @param covariates Covariate data (not used - expression already adjusted)
#' @param brain_region Brain region name
#' @param CONFIG Configuration parameters
#' @return List with model results
train_elastic_net = function(gene_id, gene_annot, expression, genotypes, snp_annot, covariates=NULL, brain_region, CONFIG) {
  # Input validation
  if(is.null(gene_id) || is.null(expression) || is.null(genotypes) || is.null(snp_annot)) {
    stop("Missing required input data")
  }
  
  if(!gene_id %in% colnames(expression)) {
    warning(sprintf("Gene %s not found in expression data", gene_id))
    return(NULL)
  }
  
  # Get fold IDs from previous step
  fold_file = get_file_path("cv_folds", brain_region=brain_region, CONFIG=CONFIG)
  # if(!file.exists(fold_file)) {
  #   stop(sprintf("CV fold file not found: %s", fold_file))
  # }
  set.seed(070892)
  generate_fold_ids = function(n_samples, n_folds) {
    n = ceiling(n_samples / n_folds)
    fold_ids = rep(1:n_folds, n)
    sample(fold_ids[1:n_samples])
  }
  train_test_fold_ids = generate_fold_ids(nrow(genotypes), CONFIG$n_folds)
  
  if(length(train_test_fold_ids) != nrow(genotypes)) {
    stop(sprintf("Mismatch between number of samples in genotypes (%d) and CV fold IDs (%d)", 
                 nrow(genotypes), length(train_test_fold_ids)))
  }
  
  n_train_test_folds = length(unique(train_test_fold_ids))
  n_k_folds          = CONFIG$n_folds
  
  # Get adjusted expression
  stopifnot(identical(rownames(expression), rownames(genotypes)))
  adj_expression           = scale(expression[, gene_id], center=T, scale = T)
  rownames(adj_expression) = rownames(expression)
  
  # Check for NA values 
  if(any(is.na(adj_expression))) {
    warning(sprintf("NA values found in expression data for gene %s", gene_id))
    return(NULL)
  }
  
  # to get performance evaluations adapted from PrediXcan
  perf_measures <- nested_cv_elastic_net_perf(as.matrix(genotypes), (adj_expression), nrow(genotypes), train_test_fold_ids, 0.5, n_k_folds)
  R2_avg <- perf_measures$R2_avg
  rmse_avg <- perf_measures$rmse_avg
  R2_sd <- perf_measures$R2_sd
  pval_est <- perf_measures$pval_est
  rho_avg <- perf_measures$rho_avg
  rho_se <- perf_measures$rho_se
  rho_zscore <- perf_measures$rho_zscore
  rho_avg_squared <- perf_measures$rho_avg_squared
  zscore_pval <- perf_measures$zscore_pval
  
 
  ##########
  #  Fit on all data 
  ########
 
  # Run Cross-Validation
  set.seed(0212)
  cv_fold_ids <- generate_fold_ids(nrow(genotypes), n_k_folds)
  fit <- tryCatch(cv.glmnet(as.matrix(genotypes),as.vector(adj_expression), nfolds = n_k_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, keep = TRUE),
                  error = function(cond) {message('Error'); message(geterrmessage()); list()})
  if (length(fit) > 0) {
    cv_R2_folds <- rep(0, n_k_folds)
    cv_corr_folds <- rep(0, n_k_folds)
    cv_zscore_folds <- rep(0, n_k_folds)
    cv_pval_folds <- rep(0, n_k_folds)
    best_lam_ind <- which.min(fit$cvm)
    for (j in 1:n_k_folds) {
      fold_idxs <- which(cv_fold_ids == j)
      adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
      cv_R2_folds[j] <- calc_R2(adj_expression[fold_idxs], adj_expr_fold_pred)
      cv_corr_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor(adj_expr_fold_pred, adj_expression[fold_idxs]), 0)
      cv_zscore_folds[j] <- atanh(cv_corr_folds[j])*sqrt(length(adj_expression[fold_idxs]) - 3) # Fisher transformation
      cv_pval_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor.test(adj_expr_fold_pred, adj_expression[fold_idxs])$p.value, runif(1))
    }
    cv_R2_avg <- mean(cv_R2_folds)
    cv_R2_sd <- sd(cv_R2_folds)
    adj_expr_pred <- predict(fit, as.matrix(genotypes), s = 'lambda.min')
    training_R2 <- calc_R2(adj_expression, adj_expr_pred)
    cv_rho_avg <- mean(cv_corr_folds)
    cv_rho_se <- sd(cv_corr_folds)
    cv_rho_avg_squared <- cv_rho_avg**2
    # Stouffer's method for combining z scores.
    cv_zscore_est <- sum(cv_zscore_folds) / sqrt(n_k_folds)
    cv_zscore_pval <- 2*pnorm(abs(cv_zscore_est), lower.tail = FALSE)
    cv_pval_est <- pchisq(-2 * sum(log(cv_pval_folds)), 2*n_k_folds, lower.tail = F)
    if (fit$nzero[best_lam_ind] > 0) {
      # number of predictors
      n_fit=fit$nzero[best_lam_ind]
      # Discuss with Yungil: use R2 of correlation between predicted and observed expression values for all samples:
      y_all=predict(fit, as.matrix(genotypes), s = 'lambda.min')
      # correlation R2
      corr_R=ifelse(sd(y_all) != 0, cor(y_all, adj_expression),0)
      pval=cor.test(y_all, adj_expression)$p.value
      corr_R2=corr_R**2
      # adjusted correlation R2
      adj_R2=1-(1-corr_R2)*(nrow(genotypes)-1)/(nrow(genotypes)-1-n_fit)
      
      weights <- fit$glmnet.fit$beta[which(fit$glmnet.fit$beta[,best_lam_ind] != 0), best_lam_ind]
      weighted_snps <- names(fit$glmnet.fit$beta[,best_lam_ind])[which(fit$glmnet.fit$beta[,best_lam_ind] != 0)]
      # Output best betas
      bestweightlist <- weighted_snps
      bestweightinfo <- snp_annot[bestweightlist,]
      bestweightinfo = bestweightinfo[,c("Chromosome","Position","Al2","Al1","SNP")]
      
      # calculate distTSS
      if (gene_annot[gene_annot$gencodeID==gene_id,"strand"]=='-') {
        bestweightinfo[,6]=gene_annot[gene_annot$gencodeID==gene_id,"end"]-bestweightinfo[,2] }
      if (gene_annot[gene_annot$gencodeID==gene_id,"strand"]=='+') {
        bestweightinfo[,6]=bestweightinfo[,2]-gene_annot[gene_annot$gencodeID==gene_id,"start"] }
      colnames(bestweightinfo)<-c('chr', 'pos', 'refAllele', 'effectAllele','rsid','distTSS')
      
      weighttable <- data.frame(bestweightinfo, weight=weights)
      weightfile <- weighttable[,c("chr","pos","distTSS","rsid","refAllele","effectAllele","weight")]
      weightfile$gene=gene_id
      weightfile$gene_name=gene_annot[gene_annot$gencodeID==gene_id,"gene_name"]
      weightfile$gene_start=gene_annot[gene_annot$gencodeID==gene_id,"start"]
      weightfile$gene_end=gene_annot[gene_annot$gencodeID==gene_id,"end"]
      weightfile$strand=gene_annot[gene_annot$gencodeID==gene_id,"strand"]
      
      weightfile$alpha=0.5
      weightfile<-weightfile[,c('gene','gene_name','gene_start','gene_end','strand','chr','pos','distTSS','rsid','refAllele','effectAllele','weight','alpha')]
     
 
  # Create results list with consistent naming
  cv_results = list(
    n_snps_in_window        = ncol(genotypes),
    n_snps_in_model         = fit$nzero[best_lam_ind], 
    best_lambda             = fit$lambda[best_lam_ind],
    inner_R2_avg            = R2_avg,
    inner_R2_sd             = R2_sd,
    inner_pval_est          = pval_est,
    inner_rho_avg           = rho_avg,
    inner_rho_se            = rho_se,
    inner_rho_zscore        = rho_zscore,
    inner_rho_avg_squared   = rho_avg_squared,
    inner_zscore_pval       = zscore_pval,
    cv_R2_gene_avg          = cv_R2_avg,
    cv_R2_gene_sd           = cv_R2_sd,
    cv_rho_gene_avg         = cv_rho_avg,
    training_R2             = training_R2,
    cv_rho_gene_se          = cv_rho_se,
    cv_rho_gene_avg_squared = cv_rho_avg_squared,
    cv_zscore_est_gene      = cv_zscore_est,
    cv_zscore_pval_gene     = cv_zscore_pval,
    cv_pval_est_gene        = cv_pval_est,
    # Performance metrics from final model
    cor_gene_all_data_pred  = corr_R,
    adj_rsq_gene            = adj_R2,
    rmse_avg                = rmse_avg,
    pval_gene               = pval,
    # Additional metrics that might be expected
    cv_adj_r2_gene_avg      = NA,
    cv_adj_r2_gene_sd       = NA,
    cv_pval_lm_gene_avg     = NA
  )
  
  result_list = list(
    cv_results         = cv_results,
    weighted_snps_info = weightfile
  )
  
  message(sprintf("Gene %s: Returning result with %d weighted SNPs", gene_id, nrow(weightfile)))
  
  return(result_list)
  }else { 
    message(sprintf("Gene %s: there are no snps selected regarding the gene", gene_id))
    return(NULL)
  }
  
  }else{
    message(sprintf("Gene %s: there is no Model", gene_id))
    return(NULL)
  }
}

