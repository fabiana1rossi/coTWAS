#####################################################################
# Module Training Core Functions
#####################################################################
# This file contains the core functions for:
# 1. Data preparation and processing
# 2. Model training and cross-validation
# 3. Performance evaluation
#####################################################################


#' Prepare feature matrix for training
#' @param genotypes Genotype data
#' @param covariates Covariate data
#' @param significant_snps Vector of significant SNPs
#' @return Feature matrix
prepare_genotype_matrix = function(genotypes, covariates, significant_snps) {
    # Subset genotypes for significant coeQTLs
    geno_subset = genotypes[, colnames(genotypes) %in% significant_snps, drop = FALSE]
    
    if(!is.null(covariates)) {
        # Combine genotypes and covariates
        feature_matrix = cbind(
            geno_subset,
            covariates
        )
    } else {
        feature_matrix = geno_subset
    }
    
    return(as.matrix(feature_matrix))
}

#' Get cis genotypes for a gene
#' @param gt_df Genotype data frame
#' @param snp_annot SNP annotation data
#' @param coords Gene coordinates [start, end]
#' @param cis_window Window size for cis-SNPs
#' @return Matrix of cis genotypes
get_cis_genotype = function(gt_df, snp_annot, coords, cis_window) {
    snp_info = snp_annot %>% dplyr::filter(Chromosome == coords[1] & Position >= (coords[2] - cis_window) & Position <= (coords[3] + cis_window))
    
    if (nrow(snp_info) == 0)
        return(NA)
    
    if (TRUE %in% (snp_info$SNP %in% colnames(gt_df))) {
        cis_gt = gt_df[, colnames(gt_df) %in% snp_info$SNP, drop = FALSE]
    } else {
        return(NA)
    }
    
    return(colnames(cis_gt))
}

#' Get gene type from annotation
get_gene_type = function(gene_annot, gene) {
    filter(gene_annot, gencodeID == gene)$gene_type
}

# Get gene name from annotation
get_gene_name = function(gene_annot, gene) {
    filter(gene_annot, gencodeID == gene)$gene_name
}

#' Get gene coordinates from annotation
get_gene_coords = function(gene_annot, gene) {
    row = gene_annot[which(gene_annot$gencodeID == gene),]
    c(row$chr, row$start, row$end)
}



#' Train elastic net model for a module
#' @param brain.region Brain region name
#' @param x Genotype matrix
#' @param y Expression matrix
#' @param gene Gene name
#' @param covs Covariates
#' @param cv_fold_id CV fold identifier
#' @return List of performance measures and model
train_elastic_net_module = function(brain.region, network_name, module_name, data_file_path,significant_coeqtls, x, y, gene, gene_annot, snp_annot,covs, cv_fold_id, cis_window = 1e6) {
    
    # load fold IDs from coEQTL previous step
    fold_file = file.path(data_file_path, 
                          paste0(cv_fold_id, '_cv_4fold_', brain.region, '_ids.RData'))
    train_test_fold_ids = get(load(fold_file))  


    ## Filter genotype matrix for selected eqtls
    x = x[, colnames(x) %in% significant_coeqtls, drop = FALSE]
    if(nrow(x) == 0){return(print('No significant SNPs'))}
  
    ## Match samples expression, genotype and covariates
    samples = intersect(rownames(y), rownames(x))
    y       = y[which(rownames(y) %in% samples), ,drop = FALSE]
    x       = x[which(rownames(x) %in% samples), colSums(is.na(x)) == 0, drop = FALSE]
    covs    = covs[which(rownames(covs) %in% samples), ]
   
    y        = y[match(rownames(x), rownames(y)),]
    covs     = covs[match(rownames(x), rownames(covs)),]

    stopifnot(identical(rownames(y), rownames(x)));stopifnot(identical(rownames(y), rownames(covs)))
    
    ## Get gene type
    gene_type   = get_gene_type(gene_annot, gene)
    ## Get gene name
    gene_name   = get_gene_name(gene_annot, gene)

    # Get gene cis-SNPs (1Mbp) and remove them from genotype matrix
    gene_annot$chr = as.numeric(gene_annot$chr); gene_annot$start = as.numeric(gene_annot$start); gene_annot$end = as.numeric(gene_annot$end)
    gene_coords    = get_gene_coords(gene_annot, gene)
    cis_gt         = get_cis_genotype(x, snp_annot, gene_coords, cis_window)
    if(all(!is.na(cis_gt))) {
        x = x[, !colnames(x) %in% cis_gt]
        if(ncol(x) < 2){return(NA)}
    }
    
    # Train model using lambda_tuning
    model_results = lambda_tuning(
        brain.region = brain.region,
        samples = samples,
        n_samples = nrow(x),
        train_test_fold_ids = train_test_fold_ids,
        x            = x,
        y            = y,
        gene         = gene,
        covs         = covs,
        cv_fold_id   = cv_fold_id,
        network_name = network_name,
        module_name  = module_name,
        gene_name    = gene_name,
        gene_type    = gene_type,
        snp_annot    = snp_annot
    )
    
    return(model_results)
}

#' Calculate R-squared
#' @param y True values
#' @param y_pred Predicted values
#' @return R-squared value
calc_R2 = function(y, y_pred) {
    1 - sum((y - y_pred)^2) / sum((y - mean(y))^2)
}

#' Train elastic net model with cross-validation
#' @param X Feature matrix
#' @param Y Response variable
#' @param fold_ids Pre-defined fold assignments
#' @return Trained model and performance metrics
train_elastic_net = function(X, Y, fold_ids) {
    # Ensure inputs are matrices
    X = as.matrix(X)
    Y = as.matrix(Y)
    
    # Fit model using cv.glmnet
    fit = cv.glmnet(
        x = X,
        y = Y,
        alpha = 0.5,  
        foldid = fold_ids,
        standardize = TRUE,
        parallel = TRUE
    )
    
    # Get predictions at optimal lambda
    predictions = predict(fit, X, s = "lambda.min")
    
    # Calculate performance metrics
    metrics = list(
        rmse = sqrt(mean((Y - predictions)^2)),
        r_squared = cor(Y, predictions)^2,
        mae = mean(abs(Y - predictions))
    )
    
    return(list(
        model = fit,
        metrics = metrics,
        predictions = predictions
    ))
}


generate_fold_ids = function(n_samples, n_folds) {
  n = ceiling(n_samples / n_folds)
  fold_ids = rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}

#' Calculate cross-validation performance
#' @param model Trained model
#' @param X Feature matrix
#' @param Y Response variable
#' @param fold_ids Pre-defined fold assignments
#' @return CV performance metrics
calculate_cv_performance = function(model, X, Y, fold_ids) {
    unique_folds = unique(fold_ids)
    cv_predictions = numeric(length(Y))
    
    for(fold in unique_folds) {
        test_idx = which(fold_ids == fold)
        cv_predictions[test_idx] = predict(model, 
                                          newx = X[test_idx,], 
                                          s = "lambda.min")
    }
    
    # Calculate metrics
    performance = list(
        rmse_cv = sqrt(mean((Y - cv_predictions)^2)),
        r_squared_cv = cor(Y, cv_predictions)^2,
        mae_cv = mean(abs(Y - cv_predictions))
    )
    
    return(performance)
}

#' Save model results
#' @param model_results Model and performance results
#' @param output_path Output file path
#' @return NULL
save_model_results = function(model_results, output_path) {
    # Create results summary
    results_summary = list(
        model = model_results$model,
        performance = list(
            training = model_results$metrics,
            cv = model_results$model$cvm  # Use cv.glmnet's built-in CV metrics
        ),
        lambda_min = model_results$model$lambda.min,
        lambda_1se = model_results$model$lambda.1se,
        feature_coefficients = coef(model_results$model, s = "lambda.min")
    )
    
    # Save results
    save(results_summary,
         file = output_path,
         compress = TRUE)
}

#' Generate model diagnostics
#' @param model_results Model and performance results
#' @param output_dir Output directory
#' @return NULL
generate_model_diagnostics = function(model_results, output_dir) {
    pdf(file.path(output_dir, "model_diagnostics.pdf"))
    
    # Plot 1: Predicted vs Actual
    plot(model_results$predictions,
         model_results$model$y,
         main = "Predicted vs Actual",
         xlab = "Predicted",
         ylab = "Actual")
    abline(0, 1, col = "red")
    
    # Plot 2: Cross-validation plot
    plot(model_results$model)
    
    # Plot 3: Coefficient path
    plot(model_results$model$glmnet.fit, xvar = "lambda")
    abline(v = log(model_results$model$lambda.min), col = "red")
    abline(v = log(model_results$model$lambda.1se), col = "blue")
    
    dev.off()
}

#' Validate model assumptions
#' @param model_results Model and performance results
#' @return List of validation results
validate_model_assumptions = function(model_results) {
    residuals = model_results$model$y - model_results$predictions
    
    # Perform validation tests
    validation_results = list(
        normality = shapiro.test(residuals),
        residual_stats = list(
            mean = mean(residuals),
            sd = sd(residuals),
            skewness = moments::skewness(residuals),
            kurtosis = moments::kurtosis(residuals)
        )
    )
    
    return(validation_results)
}

#' Get matched modules for a brain region
#' @param modules Module list
#' @param expression_data Expression data
#' @param brain_region Brain region name
#' @return Matched module data
get_matched_modules = function(modules, expression_data, brain_region) {
    # Implementation of module matching logic
    matched = modules %>%
        tibble::enframe() %>%
        tidyr::unnest_longer(value) %>%
        filter(value %in% colnames(expression_data))
    
    matched
}

#' Train models for each network
#' @param matched_modules Matched module data
#' @param expression_data Expression data
#' @param brain_region Brain region name
#' @param output_dir Output directory
#' @param n_cores Number of cores for parallel processing
#' @return Network training results
train_network_models = function(matched_modules, expression_data, brain_region, 
                                output_dir, n_cores = 1) {
    # Create output directories
    model_dir = file.path(output_dir, "models", brain_region)
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Train models in parallel
    network_results = mclapply(unique(matched_modules$name), function(network) {
        message("Training network: ", network)
        
        # Get network modules
        network_modules = matched_modules %>%
            filter(name == network)
        
        # Train elastic net model
        model = train_elastic_net(
            modules = network_modules,
            expression_data = expression_data,
            output_dir = model_dir
        )
        
        # Save model
        save(model, 
             file = file.path(model_dir, paste0(network, "_model.RData")),
             compress = TRUE)
        
        model
    }, mc.cores = n_cores)
    
    names(network_results) = unique(matched_modules$name)
    network_results
}

#' Create model databases
#' @param network_results Network training results
#' @param brain_region Brain region name
#' @param output_dir Output directory
#' @return List of model databases
create_model_databases = function(network_results, brain_region, output_dir) {
    # Implementation of database creation
    db_dir = file.path(output_dir, "databases", brain_region)
    dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Create databases for each network
    model_dbs = lapply(names(network_results), function(network) {
        create_network_database(
            model = network_results[[network]],
            network = network,
            output_dir = db_dir
        )
    })
    
    names(model_dbs) = names(network_results)
    model_dbs
}

#' Evaluate prediction performance
#' @param model_dbs Model databases
#' @param expression_data Expression data
#' @param brain_region Brain region name
#' @param output_dir Output directory
#' @return Performance metrics
evaluate_prediction_performance = function(model_dbs, expression_data, 
                                          brain_region, output_dir) {
    # Implementation of performance evaluation
    results_dir = file.path(output_dir, "results", brain_region)
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Calculate performance metrics
    performance = lapply(names(model_dbs), function(network) {
        calculate_performance_metrics(
            model_db = model_dbs[[network]],
            expression_data = expression_data,
            network = network
        )
    })
    
    names(performance) = names(model_dbs)
    
    # Save performance results
    save(performance,
         file = file.path(results_dir, "performance_metrics.RData"),
         compress = TRUE)
    
    performance
}

#' Create network database
#' @param model Trained model
#' @param network Network name
#' @param output_dir Output directory
#' @return Database object
create_network_database = function(model, network, output_dir) {
    # Implementation of database creation
}

#' Calculate performance metrics
#' @param model_db Model database
#' @param expression_data Expression data
#' @param network Network name
#' @return Performance metrics
calculate_performance_metrics = function(model_db, expression_data, network) {
    # Implementation of performance calculation
}
#' Test module eigengene with lambda tuning
#' @param brain.region Brain region name
#' @param samples Sample IDs
#' @param n_samples Number of samples
#' @param train_test_fold_ids CV fold assignments
#' @param x Feature matrix
#' @param y Expression matrix
#' @param gene Gene name
#' @param covs Covariates
#' @param cv_fold_id CV fold identifier
#' @param network_name Network name
#' @param module_name Module name
#' @param gene_name Gene name 
#' @param gene_type Gene type
#' @param snp_annot SNP annotation data
#' @return List of performance measures and model
lambda_tuning = function(brain.region, samples, n_samples, train_test_fold_ids, x, y,
                         gene, covs, cv_fold_id, network_name, module_name,
                         gene_name, gene_type, snp_annot) {
    
    # Generate lambda sequence
    lambda_max = max(abs(t(as.matrix(x)) %*% as.matrix(y)), na.rm = TRUE) / nrow(x)
    lambda_min = lambda_max * 1e-02
    lambda_seq = exp(seq(log(lambda_max), log(lambda_min), length.out = 100))
    
    # Cross validation loop
    cv_results = lapply(1:length(unique(train_test_fold_ids)), function(test_fold) {
        # Split data
        train_idxs = which(train_test_fold_ids != test_fold)
        test_idxs  = which(train_test_fold_ids == test_fold)
        
        x_train    = x[rownames(x) %in% samples[train_idxs], , drop = FALSE]
        y_train    = y[rownames(y) %in% rownames(x_train), , drop = FALSE]
        x_test     = x[rownames(x) %in% samples[test_idxs], , drop = FALSE]
        y_test     = y[rownames(y) %in% rownames(x_test), , drop = FALSE]
        
        # Calculate principal components
        train_pca    = prcomp(as.matrix(y_train), scale = TRUE, center = TRUE)
        loadings_ref = train_pca$rotation[,1]
        y_train_pc1  = train_pca$x[,1]
        
        # Project test data onto training PCs
        loadings_test = loadings_ref[colnames(y_test)]
        y_test_pc1    = as.matrix(y_test) %*% loadings_test
        
        # Fit model
        set.seed(123)
        ## Inner loop 
        cv_inner_folds = generate_fold_ids(length(y_train_pc1), 4)
        fit            = cv.glmnet(as.matrix(x_train), y_train_pc1, 
                          alpha = 0.5, foldid = cv_inner_folds,
                          lambda = lambda_seq)
        
        # Make predictions
        pred = predict(fit, newx = as.matrix(x_test), s = 'lambda.min')
        
        # Calculate performance metrics
        list(
          fit = fit,
          
          # === Primary: Performance on PC1 (the actual model target) ===
          r2_pc1          = calc_R2(y_test_pc1, pred),
          correlation_pc1 = cor(y_test_pc1, pred),
          pvalue_cor_pc1  = cor.test(y_test_pc1, pred)$p.value,
          adj_r2_pc1      = summary(lm(y_test_pc1 ~ pred))$adj.r.squared,
          pvalue_pc1      = summary(lm(y_test_pc1 ~ pred))$coefficients[2,4],
          zscore_pc1      = atanh(cor(y_test_pc1, pred)) * sqrt(length(pred) - 3),
         
          # === Secondary: Performance on the gene (optional check) ===
          r2_gene          = calc_R2(y_test[, gene], pred),
          correlation_gene = cor(y_test[, gene], pred),
          pvalue_cor_gene  = cor.test(y_test[, gene], pred)$p.value,
          adj_r2_gene      = summary(lm(y_test[, gene] ~ pred))$adj.r.squared,
          pvalue_gene      = summary(lm(y_test[, gene] ~ pred))$coefficients[2,4],
          zscore_gene      = atanh(cor(y_test[, gene], pred)) * sqrt(length(pred) - 3)
        )
    })
    
    # Find best lambda across folds
    cv_errors   = sapply(cv_results, function(r) r$fit$cvm)
    best_lambda = lambda_seq[which.min(colMeans(cv_errors))]
    
    # Fit final model
    y_pc1 = prcomp(as.matrix(y), scale = TRUE, center = TRUE)$x[,1]
    
    # Flip the sign of PC1 if it's negatively correlated with the gene expression
    if(cor(y[,gene], y_pc1) < 0) {
        sign    = TRUE
        y_pc1   = -y_pc1
    }else{sign = FALSE}
    
    final_model = cv.glmnet(as.matrix(x), y_pc1, alpha = 0.5,
                            lambda = c(best_lambda-1e-8, best_lambda))
    
    # Calculate final predictions and performance (as EpiXcan)
    predictions = predict(final_model, as.matrix(x), s = 'lambda.min')
    
    # Calculate additional test statistics
    n_train_test_folds      = length(unique(train_test_fold_ids))
    cv_zscores_gene         = sapply(cv_results, function(r) r$zscore_gene)
    cv_pvalues_gene         = sapply(cv_results, function(r) r$pvalue_gene)
    
    cv_zscores_pc1         = sapply(cv_results, function(r) r$zscore_pc1)
    cv_pvalues_pc1         = sapply(cv_results, function(r) r$pvalue_pc1)
    
    # Compile performance metrics
    performance = list(
      
        # Primary (Performance on PC1)
        cv_R2_pc1_avg          = mean(sapply(cv_results, function(r) r$r2_pc1)),
        cv_R2_pc1_sd           = sd(sapply(cv_results, function(r) r$r2_pc1)),
        cv_rho_pc1_avg         = mean(sapply(cv_results, function(r) r$correlation_pc1)),
        cv_rho_pc1_se          = sd(sapply(cv_results, function(r) r$correlation_pc1)),
        cv_rho_pc1_avg_squared = mean(sapply(cv_results, function(r) r$correlation_pc1))^2,
        cv_zscore_est_pc1      = sum(cv_zscores_pc1) / sqrt(n_train_test_folds),
        cv_zscore_pval_pc1     = 2 * pnorm(abs(sum(cv_zscores_pc1) / sqrt(n_train_test_folds)), 
                                            lower.tail = FALSE),
        cv_pval_est_pc1            = pchisq(-2 * sum(log(cv_pvalues_pc1)), 
                                         df                      = 2*n_train_test_folds, lower.tail = FALSE),
        
        # Secondary (performance on gene)
        cv_R2_gene_avg          = mean(sapply(cv_results, function(r) r$r2_gene)),
        cv_R2_gene_sd           = sd(sapply(cv_results, function(r) r$r2_gene)),
        cv_rho_gene_avg         = mean(sapply(cv_results, function(r) r$correlation_gene)),
        cv_rho_gene_se          = sd(sapply(cv_results, function(r) r$correlation_gene)),
        cv_rho_gene_avg_squared = mean(sapply(cv_results, function(r) r$correlation_gene))^2,
        cv_zscore_est_gene      = sum(cv_zscores_gene) / sqrt(n_train_test_folds),
        cv_zscore_pval_gene     = 2 * pnorm(abs(sum(cv_zscores_gene) / sqrt(n_train_test_folds)), 
                                  lower.tail = FALSE),
        cv_pval_est_gene        = pchisq(-2 * sum(log(cv_pvalues_gene)), 
        df                      = 2*n_train_test_folds, lower.tail = FALSE),
        
        
        # Performance on all data (outside inner cv) 
        cor_pc1_all_data_pred       = cor(y_pc1, predictions),
        cor_gene_all_data_pred      = cor(y[,gene], predictions),
        
        adj_rsq_gene_all_data       = summary(lm(y[,gene] ~ predictions))$adj.r.squared,
        adj_rsq_pc1_all_data        = summary(lm(y_pc1 ~ predictions))$adj.r.squared,
        pval_gene_all_data          = tryCatch(summary(lm(y[,gene] ~ predictions))$coef[2,4], 
                                        error = function(e) NA),
        pval_pc1_all_data           = tryCatch(summary(lm(y_pc1 ~ predictions))$coef[2,4], 
                                        error = function(e) NA),
        sign = sign
    )
    
    # Extract weights if model converged
    if (length(final_model) > 0) {
        # Get weights from final model
        weights = data.frame(snp = rownames(coef(final_model, s = "lambda.min"))[-1], weights = as.numeric(coef(final_model, s = "lambda.min"))[-1])  # Exclude intercept
        weights = weights[weights$weights != 0, ]  # Keep only non-zero weights
        weighted_snps = weights$snp
        
        # Create weighted SNPs info dataframe
        weighted_snps_info = snp_annot %>% 
            filter(SNP %in% weighted_snps) %>% 
            dplyr::select(SNP, varID, Al2, Al1) %>% dplyr::rename(snp=SNP, ref_vcf = Al2, alt_vcf = Al1) %>%
            mutate(gene = gene) %>%
            merge(., weights, by = "snp") %>%
            dplyr::select(gene, snp, varID, ref_vcf, alt_vcf, weights)
        
        # Create model summary
        model_summary = c(
            gene, gene_name, gene_type,
            network_name,
            module_name,
            0.5,  # alpha value
            ncol(x),
            length(unique(weighted_snps)),
            final_model$lambda.min,
            performance$cv_R2_pc1_avg,
            performance$cv_R2_pc1_sd,  
            performance$cv_rho_pc1_avg,
            performance$cv_rho_pc1_se,
            performance$cv_rho_pc1_avg_squared,
            performance$cv_zscore_est_pc1,
            performance$cv_zscore_pval_pc1,
            performance$cv_pval_est_pc1,
            performance$cv_R2_gene_avg,
            performance$cv_R2_gene_sd,
            performance$cv_rho_gene_avg,
            performance$cv_rho_gene_se,
            performance$cv_rho_gene_avg_squared,
            performance$cv_zscore_est_gene,
            performance$cv_zscore_pval_gene,
            performance$cv_pval_est_gene,
            performance$cor_pc1_all_data_pred,
            performance$cor_gene_all_data_pred,
            performance$adj_rsq_gene_all_data,
            performance$adj_rsq_pc1_all_data,
            performance$pval_gene_all_data,
            performance$pval_pc1_all_data,
            performance$sig
        )
        
       
        return(list(
            perf.measures = performance,
            model = final_model,
            weighted_snps_info = weighted_snps_info,
            model_summary = model_summary
        ))
    } else {
        # Create empty weighted SNPs info for failed model
        weighted_snps_info = as.data.frame(matrix(NA, nrow=1, ncol=6))
        colnames(weighted_snps_info) = c("gene", "snp", "varID", "ref_vcf", "alt_vcf", "weights")
        
        # Create model summary with default values for failed model
        model_summary = c(
            gene, gene_name, gene_type,
            network_name, module_name,
            0.5,  # alpha value
            ncol(x),
            0,    # ind
            0,    # lambda
            NA,
            NA,   
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,   
            NA,
            NA, NA, NA, NA,  
            NA,
            sign
        )
        
        return(list(
            perf.measures = performance,
            model = final_model,
            weighted_snps_info = weighted_snps_info,
            model_summary = model_summary
        ))
    }
}
