#################################################
## Cross-Validation Helper Functions
#################################################

#' Generate stratified fold IDs based on diagnosis
#' @param n_folds Number of folds
#' @param covariates Covariate data with Dx column
#' @return Vector of stratified fold assignments
generate_fold_ids_stratified = function(n_folds, covariates) {
  set.seed(070892)
  covariates$Dx = as.character(covariates$Dx)
  unique_dx = unique(covariates$Dx)
  
  if(length(table(covariates$Dx)) == 1) {
    return(generate_fold_ids_notstratified(n_folds, covariates))
  }
  
  fold_assignments = map(unique_dx, function(dx_value) {
    dx_data = covariates %>% filter(Dx == dx_value)
    sample(1:n_folds, size = nrow(dx_data), replace = TRUE)
  })
  
  do.call(c, fold_assignments)
}

#################################################
## Statistical Analysis Functions
#################################################

#' Perform RLM analysis on genotype and PC1 data
#' @param x_train Genotype data
#' @param y_train_pc1 Expression PC1 data
#' @return List of RLM results
perform_rlm_analysis = function(x_train, y_train_pc1, n_cores=1) { ## increase n_cores to run in parallel 
  fit_rlm = function(G, y) {
    rlm(y ~ G, psi = "psi.bisquare")
  }
  
  QTL.rlm = mclapply(colnames(x_train), function(snp) {
    if(length(unique(x_train[,snp])) == 1) return(NA)
    fit_rlm(G = x_train[,snp], y = y_train_pc1[,1])
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  
  setNames(QTL.rlm[!is.na(QTL.rlm)], colnames(x_train)[!is.na(QTL.rlm)])
}

#' Process RLM results to get statistics
#' @param QTL.rlm List of RLM results
#' @return Data frame of statistics
process_rlm_results = function(QTL.rlm, n_cores = 1) { ## increase number of cores to run in parallel
  # Extract coefficients and p-values
  coefficients = pbsapply(QTL.rlm, function(l) summary(l)$coeff[2,], cl = 1)
  p_values = unlist(mclapply(QTL.rlm, function(m) f.robftest(m, var = "G")$p.value, 
                             mc.cores = n_cores, mc.preschedule = FALSE))
  
  # Combine results
  results = data.frame(
    t(coefficients),
    p_value = p_values,
    stringsAsFactors = FALSE
  )
  colnames(results) = c("beta", "st_err", "t_value", "p_value")
  
  # Add ranks
  results = results[order(results$p_value),]
  results$rank = rank(results$p_value)
  
  results
}

#################################################
## LD Pruning Functions
#################################################

#' Calculate correlation and prune SNPs
#' @param genotypes Genotype data
#' @param flank Flanking region size
#' @param association Association data
#' @param snp SNP index
#' @param rsquared R-squared threshold
#' @return List of pruning results
SNPr2 = function(genotypes, flank, association, snp, rsquared) {
  # Find SNPs in the flanking region
  index = which(association[,2] == association[snp,2] &
                  abs(association[,3] - association[snp,3]) <= flank)
  
  # Calculate LD
  marker = row.names(association)[index]
  window = genotypes[, colnames(genotypes) %in% marker, drop = FALSE]
  r2 = cor(window, genotypes[,association[snp,1]], use = "pairwise.complete.obs")^2
  
  # Identify SNPs to prune
  snp_r2 = row.names(subset(r2, r2[,1] >= rsquared))
  snp_r2 = snp_r2[snp_r2 != association[snp,1]]
  
  list(
    cut = subset(r2, r2[,1] >= rsquared),
    new.genotypes = genotypes[, !colnames(genotypes) %in% snp_r2],
    new.association = association[!association[,1] %in% snp_r2,],
    snp = snp_r2
  )
}


prune_coeqtls = function(annot ,coeqtls,x, f, rsq){
  message("Performing ranked pruning")
  
  annot = annot[which(rownames(annot) %in% rownames(coeqtls)),]
  annot = annot[match(rownames(coeqtls),rownames(annot)),]
  
  assoc           = cbind(annot[, c("SNP", "Chromosome", "Position")], coeqtls[,c("new_rank", "new_pval")])
  colnames(assoc) = c("Marker", "Chromosome", "Position", "new_rank", "new_pval")
  assoc           = assoc[order(assoc$new_rank),]
  
  for(i in 1:nrow(assoc)) {
    if (i <= nrow(assoc)){
      pruning = SNPr2(genotypes    = x,
                      flank       = f,
                      association = assoc,
                      snp         = i,
                      rsquared    = rsq)
      x     = pruning$new.genotypes
      assoc = pruning$new.association
    } else { break }
  }
  
  ## Assign new rank 
  assoc$rank = rank(assoc[,"new_pval"])
  
  return(assoc)
  
}

#################################################
## Results Processing Functions
#################################################

#' Merge and multiply ranks across folds
#' @param dataframes List of dataframes with rank information
#' @return Merged dataframe with new ranks
merge_and_multiply_ranks = function(dataframes) {
  merged_df = do.call("cbind", dataframes)
  rownames(merged_df) = merged_df$SNP
  
  rank_cols = merged_df[, grepl('rank', colnames(merged_df))]
  pval_cols = merged_df[, grepl('p_value', colnames(merged_df))]
  
  results = data.frame(
    new_rank = rowProds(as.matrix(rank_cols), useNames = TRUE),
    new_pval = rowProds(as.matrix(pval_cols), useNames = TRUE)
  )
  
  # Extract beta and st_err columns for meta-analysis
  beta_cols <- merged_df[, grepl('^beta(\\.|$)', colnames(merged_df)), drop = FALSE]
  se_cols <- merged_df[, grepl('^st_err(\\.|$)', colnames(merged_df)), drop = FALSE]
  
  # Apply meta-analysis row-wise
  meta_results <- t(apply(seq_len(nrow(beta_cols)), 1, function(i) {
    betas <- as.numeric(beta_cols[i, ])
    ses <- as.numeric(se_cols[i, ])
    meta_effect(betas, ses)
  }))
  colnames(meta_results) <- c("meta_beta", "meta_se", "meta_p")
  
  # Combine meta-analysis results with merged_df
  merged_df <- cbind(merged_df, meta_results)
  
  # Sort by meta_p
  merged_df <- merged_df[order(merged_df$meta_p), ]
  
  results
}


#' Process results from all folds
#' @param coeQTLs_folds List of fold results
#' @param common_coeQTLs Vector of common SNPs
#' @return List of processed results
process_fold_results = function(coeQTLs_folds, common_coeQTLs) {
  map(coeQTLs_folds, function(fold_result) {
    result = fold_result[rownames(fold_result) %in% common_coeQTLs,]
    result = result[match(common_coeQTLs, rownames(result)),]
    result$SNP = rownames(result)
    result
  })
}

#################################################
## Main Analysis Functions
#################################################

#' Perform cross-validation for coeQTL detection
#' @param brain.region Brain region being analyzed
#' @param module_data Module data containing genotype and eigenvector
#' @param n_train_test_folds Number of CV folds
#' @param covs Covariate data
#' @param cv_fold_id Identifier for CV folds
#' @return List of coeQTL results
perform_coeqtl_cv = function(brain.region, data_file_path, module_data, n_train_test_folds,expr_module,  
                             covs, cv_fold_id) {
  
  # Remove all SNPs containing any NAs ---> missing values not supported by the Elastic net in the training step 
  genotype = module_data$geno.mod
  genotype = genotype[, colSums(is.na(genotype)) == 0, drop=FALSE]
  
  ## match genotype, expression and covariates data samples 
  samples  = intersect(rownames(expr_module),rownames(genotype))
  samples  = intersect(rownames(covs),samples)
  
  expr_module = expr_module[which(rownames(expr_module) %in% samples),]
  genotype    = genotype[which(rownames(genotype) %in% samples), colSums(is.na(genotype)) == 0]
  expr_module = expr_module[match(rownames(genotype), rownames(expr_module)),]
  covs        = covs[which(rownames(covs) %in% samples),]
  covs        = covs[match(rownames(genotype), rownames(covs)),]
  
  stopifnot(identical(rownames(covs),rownames(expr_module))); stopifnot(identical(rownames(covs),rownames(genotype)))
  
  # Generate or load fold IDs
  fold_file = file.path(data_file_path, 
                        paste0(cv_fold_id, '_cv_4fold_', brain.region, '_ids.RData'))
  
  train_test_fold_ids = if(!file.exists(fold_file)) {
    ids = generate_fold_ids_stratified(n_train_test_folds, covariates = covs)
    save(ids, file = fold_file)
    ids
  } else {
    get(load(fold_file))
  }
  
  
  # Perform cross-validation
  coeQTLs_folds = map(1:n_train_test_folds, function(test_fold) {
    # Split data
    train_idxs  = which(train_test_fold_ids != test_fold)
    x_train     = genotype[train_idxs, ]
    y_train     = expr_module[train_idxs, , drop = FALSE]
    # Compute pc1 fold.
    y_train_pc1 = as.data.frame(prcomp(y_train, scale=T, center=T)$x[,1])
    
    stopifnot(identical(rownames(x_train),rownames(y_train_pc1)))
    
    # Analyze
    QTL.rlm = perform_rlm_analysis(x_train, y_train_pc1)
    process_rlm_results(QTL.rlm)
  })
  
  # Process results: Get eqtls replicating in all folds
  common_coeQTLs = Reduce(intersect, map(coeQTLs_folds, ~ rownames(.x)))
  coeQTLs        = process_fold_results(coeQTLs_folds, common_coeQTLs)
  
  list(
    coeQTLs = coeQTLs,
    coeQTLs_newrank = merge_and_multiply_ranks(coeQTLs)
  )
}

#' Process a single module for coeQTL analysis
#' @param network_file Path to module/network RData file
#' @param covs Covariate data
#' @param n_folds Number of CV folds
#' @return Updated module results
process_module = function(network_file, data_dir, brain_region, covs, expr,  n_folds, flank_size, r_squared_threshold) {
  # Load network data
  res_module = get(load(network_file))
  # Check if computation is already performed
  if(length(res_module)==7){return(NULL)}
  
  # Process each module
  module_genes = res_module[["genes"]]
  message("Processing file: ", basename(network_file))
  
  ## Subset expression with module genes
  expr_module  = expr[,which(colnames(expr) %in% module_genes),drop=F]
  
  # Perform coeQTL analysis
  coeqtl_results   = perform_coeqtl_cv(
    brain.region       = brain_region,
    data_file_path     = data_dir,
    module_data        = res_module,
    n_train_test_folds = n_folds,
    covs               = covs,
    expr_module        = expr_module, 
    cv_fold_id         = "MODULE"
  )
  
  # Prune results if needed
  if(!is.null(coeqtl_results$coeQTLs_newrank)) {
    pruned_results   = prune_coeqtls(
      annot               = res_module$map.mod,
      coeqtls             = coeqtl_results$coeQTLs_newrank,
      x                   = res_module$geno.mod,
      f                   = flank_size,
      rsq                 = r_squared_threshold
    )
    
    res_module$coeqtls_pruned = pruned_results
    res_module$coeqtl_folds   = coeqtl_results
  }
  
  ## Update file
  save(res_module, file = network_file)
  rm(res_module); gc()
  
  return(NULL)

}

#' Run the complete coeQTL pipeline
#' @param network_res_dir Network results directory
#' @param brain_regions Brain regions to process
#' @param n_cores Number of cores for parallel processing
#' @param n_folds Number of CV folds
#' @param covs Covariate data
#' @return NULL
run_coeqtl_pipeline = function(expr_data_covs, network_res_dir,data_dir, brain_regions,
                               n_cores, n_folds, flank_size, r_squared_threshold) {
  
  for(brain_region in brain_regions) {
    message("Processing brain region: ", brain_region)
    
    ## Get covs
    training_data = get(load(expr_data_covs))
    covs          = as.data.frame(training_data[[brain_region]]$colData)
    expr          = as.data.frame(training_data[[brain_region]]$assays$expression)
    
    covs          = covs[match(rownames(expr),rownames(covs)),]
    stopifnot(identical(rownames(covs),rownames(expr)))
    
    ## Get network files
    network_files = list.files(
      network_res_dir,
      pattern = paste0(brain_region, "_network"),
      full.names = TRUE
    )
    
    
    # Process modules
    all_networks_results = mclapply(network_files, function(network_file) { 
      
      network_name = gsub(".RData","",strsplit2(basename(network_file),"network_")[,2]) 
      # Getting coeQTL and pruning
      process_module(network_file,data_dir, brain_region, covs,expr, n_folds, flank_size, r_squared_threshold)
      
      # names(results) = network_name
      # # Save results
      # save(results, 
      #      file = network_file)
      
    }, mc.cores = n_cores, mc.preschedule = FALSE)
    
    
    
    gc(verbose = FALSE)
  }
}

# Function for fixed-effect inverse variance meta-analysis
meta_effect <- function(betas, ses) {
  weights <- 1 / (ses^2)
  beta_meta <- sum(weights * betas) / sum(weights)
  se_meta <- sqrt(1 / sum(weights))
  z_meta <- beta_meta / se_meta
  p_meta <- 2 * pnorm(-abs(z_meta))
  return(c(beta_meta, se_meta, p_meta))
}
