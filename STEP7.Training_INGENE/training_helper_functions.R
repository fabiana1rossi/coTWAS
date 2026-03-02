#####################################################################
# INGENE Training Helper Functions
#####################################################################
# This file contains helper functions for the INGENE training pipeline.
# Key components:
# - Cross-validation and model evaluation functions
# - Elastic net model training functions with comprehensive metrics
# - Performance metric calculations matching CIS/EpiXcan standards
# - File handling utilities
# - Model summary and weights columns matching CIS/EpiXcan format
# - Uses imputed expression of co-expression partners as predictors (not SNPs)
#####################################################################

#' Calculate R-squared value between predicted and actual values
#' @param y True/observed values
#' @param y_pred Predicted values
#' @return R-squared value (proportion of variance explained)
calc_R2 = function(y, y_pred) {
    tss = sum(y**2)  # Total sum of squares
    rss = sum((y - y_pred)**2)  # Residual sum of squares
    1 - rss/tss  # R-squared formula
}

#' Calculate correlation between predicted and actual values
#' @param y True/observed values
#' @param y_pred Predicted values
#' @return Correlation coefficient
calc_corr = function(y, y_pred) {
    sum(y*y_pred) / (sqrt(sum(y**2)) * sqrt(sum(y_pred**2)))
}

#' Generate cross-validation fold assignments
#' @param n_samples Number of samples to split into folds
#' @param n_folds Number of folds (default: 4)
#' @return Vector of randomized fold assignments
generate_fold_ids = function(n_samples, n_folds=4) {
    n = ceiling(n_samples / n_folds)
    fold_ids = rep(1:n_folds, n)
    sample(fold_ids[1:n_samples])
}

#' Adjust expression data for covariates
#' @param expression Expression data matrix
#' @param cov_df Covariates data frame
#' @return Covariate-adjusted expression data
adjust_for_covariates = function(expression, cov_df) {
    if(!is.null(cov_df)){
        print("Correcting expression for covariates")
        cov_df = cov_df[,!grepl("BrNum|RNum|Region|Race",colnames(cov_df))]
        samples = rownames(cov_df)
        
        # Convert numeric columns
        numeric_cols = c("RIN", "mitoRate", "Age", "rRNA_rate", "overallMapRate", "totalAssignedGene")
        for(col in numeric_cols) {
            if(col %in% colnames(cov_df)) {
                cov_df[[col]] = as.numeric(as.character(cov_df[[col]]))
            }
        }
        
        # Convert factor columns
        factor_cols = c("Protocol", "Dataset")
        for(col in factor_cols) {
            if(col %in% colnames(cov_df)) {
                cov_df[[col]] = as.factor(cov_df[[col]])
            }
        }
        
        # Convert remaining columns to numeric
        if(ncol(cov_df) > 10) {
            cov_df[,11:ncol(cov_df)] = apply(cov_df[,11:ncol(cov_df)], 2, as.numeric)
        }
        
        cov_df = as.data.frame(lapply(cov_df, function(x) if (is.character(x)) as.factor(x) else x))
        rownames(cov_df) = samples
        
        stopifnot(identical(rownames(cov_df), rownames(expression)))
        
        # Adjust each gene for covariates
        expr_resid = expression
        for(i in 1:ncol(expression)) {
            combined_df = data.frame(expression_vec = expression[,i], cov_df)
            fit = lm(expression_vec ~ ., data = combined_df)
            expr_resid[,i] = residuals(fit)
        }
    } else {
        expr_resid = expression
    }
    
    # Scale the residuals
    expr_resid = scale(expr_resid, center = TRUE, scale = TRUE)
    return(expr_resid)
}

#' Get gene information from annotation
#' @param gene_annot Gene annotation data frame
#' @param gene Gene ID to look up
#' @return List with gene information
get_gene_info = function(gene_annot, gene) {
    gene_data = gene_annot[gene_annot$gencodeID == gene,]
    if(nrow(gene_data) == 0) {
        warning(sprintf("No annotation found for gene %s", gene))
        return(list(
            name = NA,
            type = NA,
            chr = NA,
            start = NA,
            end = NA,
            strand = NA
        ))
    }
    
    list(
        name = gene_data$gene_name[1],
        type = gene_data$gene_type[1],
        chr = gene_data$chr[1],
        start = gene_data$start[1],
        end = gene_data$end[1],
        strand = gene_data$strand[1]
    )
}

#' Evaluate elastic net model performance using nested cross-validation
#' @description Performs nested cross-validation to get unbiased performance estimates:
#' 1. Outer loop: Splits data into train/test folds
#' 2. Inner loop: Uses k-fold CV on training data to select lambda
#' 3. Evaluates final model on held-out test fold
#' 4. Combines results across all outer folds
#' 
#' @param x Feature matrix
#' @param y Response vector
#' @param n_samples Number of samples
#' @param train_test_fold_ids Outer CV fold assignments
#' @param alpha Elastic net mixing parameter
#' @param n_k_folds Number of inner CV folds
#' @return List of comprehensive performance metrics
nested_cv_elastic_net_perf = function(x, y, n_samples, train_test_fold_ids, alpha, n_k_folds=4) {
    n_train_test_folds = length(unique(train_test_fold_ids))
    rmse_folds = rep(0, n_train_test_folds)
    R2_folds = rep(0, n_train_test_folds)
    corr_folds = rep(0, n_train_test_folds)
    zscore_folds = rep(0, n_train_test_folds)
    pval_folds = rep(0, n_train_test_folds)
    
    # Outer-loop split into training and test set
    y = as.vector(y)
    set.seed(070892)
    
    for(test_fold in 1:n_train_test_folds) {
        train_idxs = which(train_test_fold_ids != test_fold)
        test_idxs = which(train_test_fold_ids == test_fold)
        
        x_train = x[train_idxs, ]
        y_train = y[train_idxs]
        x_test = x[test_idxs, ]
        y_test = y[test_idxs]
        
        # Inner-loop - split up training set for cross-validation to choose lambda
        cv_fold_ids = generate_fold_ids(length(y_train), n_k_folds)
        
        y_pred = tryCatch({
            # Fit model with training data
            fit = cv.glmnet(x_train, y_train, nfolds=n_k_folds, 
                           alpha=alpha, type.measure='mse', 
                           foldid=cv_fold_ids)
            # Predict test data using model that had minimal mean-squared error in cross validation
            predict(fit, x_test, s='lambda.min')
        }, error = function(cond) {
            # If the elastic-net model did not converge, predict the mean of the y_train
            rep(mean(y_train), length(y_test))
        })
        
        R2_folds[test_fold] = calc_R2(y_test, y_pred)
        
        # Get RMSE
        res = summary(lm(y_test~y_pred))
        rmse_folds[test_fold] = sqrt(1-res$r.squared)*(res$sigma)
        
        # Get p-value for correlation test between predicted y and actual y
        corr_folds[test_fold] = ifelse(sd(y_pred) != 0, cor(y_pred, y_test), 0)
        zscore_folds[test_fold] = atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3)
        pval_folds[test_fold] = ifelse(sd(y_pred) != 0, 
                                      cor.test(y_pred, y_test)$p.value, 
                                      runif(1))
    }
    
    # Aggregate results
    rmse_avg = mean(rmse_folds)
    R2_avg = mean(R2_folds)
    R2_sd = sd(R2_folds)
    rho_avg = mean(corr_folds)
    rho_se = sd(corr_folds)
    rho_avg_squared = rho_avg**2
    
    # Stouffer's method for combining z scores
    zscore_est = sum(zscore_folds) / sqrt(n_train_test_folds)
    zscore_pval = 2*pnorm(abs(zscore_est), lower.tail=FALSE)
    
    # Fisher's method for combining p-values
    pval_est = pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail=FALSE)
    
    list(
        rmse_avg=rmse_avg,
        R2_avg=R2_avg, 
        R2_sd=R2_sd, 
        pval_est=pval_est, 
        rho_avg=rho_avg, 
        rho_se=rho_se, 
        rho_zscore=zscore_est, 
        rho_avg_squared=rho_avg_squared, 
        zscore_pval=zscore_pval
    )
}

#' Train elastic net models for genes in a network
#' @description Main function for training gene expression prediction models:
#' 1. Processes each network module
#' 2. Trains models for genes in each module
#' 3. Evaluates model performance with comprehensive metrics
#' 4. Saves consolidated results with CIS/EpiXcan-compatible format
#' 
#' @param net_name Network name
#' @param brain.region Brain region name
#' @param net.data Network module definitions
#' @param train.expr Real expression data
#' @param imputed.train Imputed expression data
#' @param geneAnnot Gene annotations
#' @param covs_df Covariates data frame (optional)
#' @param output_files List of output file paths
#' @param filter Custom filter function
#' @param include.cis Whether to include cis variants (not used in INGENE)
#' @param no_cores Number of cores for parallel processing
#' @return TRUE if successful
INGENE.EnTraining = function(net_name, brain.region, net.data, train.expr,
                            imputed.train, geneAnnot, covs_df=NULL,
                            output_files, filter=NULL, include.cis=FALSE,
                            no_cores=1) {
    
    # Model parameters
    params = list(
        seed = 070891,
        alpha = 0.5,
        n_folds = 4,
        n_train_test_folds = 4
    )
    
    # Get module names for this network
    mod_names = names(net.data[[net_name]])
    
    # Initialize results storage
    all_summaries = list()
    all_weights = list()

    # Process each module
    for(mod_name in mod_names) {
        message(sprintf("Processing module: %s", mod_name))
        
        # Get genes in this module
        mod.genes = net.data[[net_name]][[mod_name]]
        
        # Subset data for module genes
        module_data = prepare_module_data(
            genes         = mod.genes,
            geneAnnot     = geneAnnot,
            imputed.train = imputed.train,
            train.expr    = train.expr,
            covs_df       = covs_df
        )
        
        if(length(module_data$mod.genes) == 0) {
            message(sprintf("No genes found for module %s, skipping", mod_name))
            next
        }
        
        # Train models for each gene with error handling
        module_results = mclapply(module_data$mod.genes, function(gene) {
            tryCatch({
                train_gene_model(
                    gene        = gene,
                    module_data = module_data,
                    net_name    = net_name,
                    mod_name    = mod_name,
                    params      = params
                )
            }, error = function(e) {
                message(sprintf("Error processing gene %s in module %s: %s", gene, mod_name, e$message))
                # Return empty results for failed genes
                list(
                    summary = create_empty_model_summary(gene, net_name, mod_name),
                    weights = data.frame()
                )
            }, warning = function(w) {
                message(sprintf("Warning processing gene %s in module %s: %s", gene, mod_name, w$message))
                # Continue with the result even if there's a warning
                tryCatch({
                    train_gene_model(
                        gene        = gene,
                        module_data = module_data,
                        net_name    = net_name,
                        mod_name    = mod_name,
                        params      = params
                    )
                }, error = function(e) {
                    message(sprintf("Error after warning for gene %s in module %s: %s", gene, mod_name, e$message))
                    list(
                        summary = create_empty_model_summary(gene, net_name, mod_name),
                        weights = data.frame()
                    )
                })
            })
        }, mc.cores = no_cores)
        
        # Combine module results with error handling
        module_summaries = tryCatch({
            # Filter out NULL results and combine summaries
            valid_summaries = lapply(module_results, function(result) {
                if(is.null(result) || is.null(result$summary) || nrow(result$summary) == 0) {
                    return(NULL)
                }
                return(result$summary)
            })
            valid_summaries = valid_summaries[!sapply(valid_summaries, is.null)]
            
            if(length(valid_summaries) > 0) {
                do.call(rbind, valid_summaries)
            } else {
                data.frame()  # Return empty data frame if no valid summaries
            }
        }, error = function(e) {
            message(sprintf("Error combining summaries for module %s: %s", mod_name, e$message))
            data.frame()  # Return empty data frame on error
        })
        
        module_weights = tryCatch({
            # Filter out NULL results and combine weights
            valid_weights = lapply(module_results, function(result) {
                if(is.null(result) || is.null(result$weights) || nrow(result$weights) == 0) {
                    return(NULL)
                }
                return(result$weights)
            })
            valid_weights = valid_weights[!sapply(valid_weights, is.null)]
            
            if(length(valid_weights) > 0) {
                do.call(rbind, valid_weights)
            } else {
                data.frame()  # Return empty data frame if no valid weights
            }
        }, error = function(e) {
            message(sprintf("Error combining weights for module %s: %s", mod_name, e$message))
            data.frame()  # Return empty data frame on error
        })
        
        all_summaries[[mod_name]] = module_summaries
        all_weights[[mod_name]] = module_weights
    }
    
    # Combine all results with error handling
    final_summaries = tryCatch({
        # Filter out empty data frames
        valid_summaries = all_summaries[sapply(all_summaries, function(x) !is.null(x) && nrow(x) > 0)]
        
        if(length(valid_summaries) > 0) {
            do.call(rbind, valid_summaries)
        } else {
            data.frame()  # Return empty data frame if no valid summaries
        }
    }, error = function(e) {
        message(sprintf("Error combining all summaries: %s", e$message))
        data.frame()  # Return empty data frame on error
    })
    
    final_weights = tryCatch({
        # Filter out empty data frames
        valid_weights = all_weights[sapply(all_weights, function(x) !is.null(x) && nrow(x) > 0)]
        
        if(length(valid_weights) > 0) {
            do.call(rbind, valid_weights)
        } else {
            data.frame()  # Return empty data frame if no valid weights
        }
    }, error = function(e) {
        message(sprintf("Error combining all weights: %s", e$message))
        data.frame()  # Return empty data frame on error
    })
    
    # Save consolidated results
    save_network_results(final_summaries, final_weights, output_files)
    
    return(TRUE)
}

#' Prepare data for a module's genes
#' @param genes Vector of gene IDs
#' @param geneAnnot Gene annotations
#' @param imputed.train Imputed expression data
#' @param train.expr Real expression data
#' @param covs_df Covariates data frame (optional)
#' @return List of prepared module data
prepare_module_data = function(genes, geneAnnot, imputed.train, train.expr, covs_df=NULL) {
    ## get genes in true df
    genes        = genes[genes %in% colnames(train.expr)]
    common.genes = intersect(colnames(imputed.train), colnames(train.expr))
    
    # Subset module genes based on genes present in both datasets
    #mod.genes = genes[which(genes %in% common.genes)]
    mod.genes = genes
    mod.genes.imputed = mod.genes[which(mod.genes %in% colnames(imputed.train))]
    
    # Adjust expression for covariates if provided
    if(!is.null(covs_df)) {
        train.expr = adjust_for_covariates(train.expr, covs_df)
    }
    
    list(
        mod.genes         = mod.genes,
        mod.genes.imputed = mod.genes.imputed,
        annotations       = geneAnnot[geneAnnot$gencodeID %in% genes,],
        imputed           = imputed.train[,colnames(imputed.train) %in% mod.genes.imputed, drop=F],
        real              = train.expr[,colnames(train.expr) %in% mod.genes, drop=F]
    )
}

#' Train model for a single gene
#' @param gene Gene ID
#' @param module_data Module data list
#' @param net_name Network name
#' @param mod_name Module name
#' @param params Model parameters
#' @return List of model results with comprehensive metrics
train_gene_model = function(gene, module_data, net_name, mod_name, params) {
    tryCatch({
        # Prepare training data
        train_data = prepare_gene_data(gene, module_data)
        
        if(is.null(train_data) || ncol(train_data$predictors) < 2) {
            return(list(
                summary = create_empty_model_summary(gene, net_name, mod_name),
                weights = data.frame()
            ))
        }
        
        # Train model with nested CV
        perf = train_with_nested_cv(train_data, params)
        
        # Create comprehensive model summary
        summary = create_comprehensive_model_summary(
            gene = gene,
            train_data = train_data,
            perf = perf,
            net_name = net_name,
            mod_name = mod_name,
            params = params
        )
        
        # Extract model weights (pass the correct weights)
        weights = extract_comprehensive_model_weights(
            predictor_weights = perf$predictor_weights,
            train_data = train_data,
            gene = gene,
            net_name = net_name,
            mod_name = mod_name,
            gene_info = get_gene_info(module_data$annotations, gene)
        )
        
        list(summary = summary, weights = weights)
        
    }, error = function(e) {
        message(sprintf("Error in train_gene_model for gene %s: %s", gene, e$message))
        # Return empty results on error
        list(
            summary = create_empty_model_summary(gene, net_name, mod_name),
            weights = data.frame()
        )
    })
}

#' Train model with nested cross-validation
#' @param train_data Training data list
#' @param params Model parameters
#' @return List of performance metrics
train_with_nested_cv = function(train_data, params) {
    set.seed(params$seed)
    
    # Generate outer CV folds
    train_test_fold_ids = generate_fold_ids(nrow(train_data$predictors), 
                                           params$n_train_test_folds)
    
    # Run nested CV
    nested_perf = nested_cv_elastic_net_perf(
        x = as.matrix(train_data$predictors),
        y = as.vector(train_data$response),
        n_samples = nrow(train_data$predictors),
        train_test_fold_ids = train_test_fold_ids,
        alpha = params$alpha,
        n_k_folds = params$n_folds
    )
    
    ##########
    #  Fit on all data 
    ########
    
    # Run Cross-Validation
    set.seed(0212)
    cv_fold_ids = generate_fold_ids(nrow(train_data$predictors), params$n_folds)
    fit = tryCatch({
        cv.glmnet(as.matrix(train_data$predictors), 
                  as.vector(train_data$response), 
                  nfolds = params$n_folds, 
                  alpha = params$alpha, 
                 # type.measure = 'mse', 
                  foldid = cv_fold_ids, 
                  keep = TRUE)
    }, error = function(cond) {
        message('Error in model training: ', geterrmessage())
        list()
    })
    
    if(length(fit) > 0) {
        cv_R2_folds = rep(0, params$n_folds)
        cv_corr_folds = rep(0, params$n_folds)
        cv_zscore_folds = rep(0, params$n_folds)
        cv_pval_folds = rep(0, params$n_folds)
        best_lam_ind = which.min(fit$cvm)
        
        for(j in 1:params$n_folds) {
            fold_idxs = which(cv_fold_ids == j)
            adj_expr_fold_pred = fit$fit.preval[fold_idxs, best_lam_ind]
            cv_R2_folds[j] = calc_R2(train_data$response[fold_idxs], adj_expr_fold_pred)
            cv_corr_folds[j] = ifelse(sd(adj_expr_fold_pred) != 0, 
                                     cor(adj_expr_fold_pred, train_data$response[fold_idxs]), 0)
            cv_zscore_folds[j] = atanh(cv_corr_folds[j]) * sqrt(length(train_data$response[fold_idxs]) - 3)
            cv_pval_folds[j] = ifelse(sd(adj_expr_fold_pred) != 0, 
                                     cor.test(adj_expr_fold_pred, train_data$response[fold_idxs])$p.value, 
                                     runif(1))
        }
        
        cv_R2_avg = mean(cv_R2_folds)
        cv_R2_sd = sd(cv_R2_folds)
        adj_expr_pred = predict(fit, as.matrix(train_data$predictors), s = 'lambda.min')
        training_R2 = calc_R2(train_data$response, adj_expr_pred)
        cv_rho_avg = mean(cv_corr_folds)
        cv_rho_se = sd(cv_corr_folds)
        cv_rho_avg_squared = cv_rho_avg**2
        
        # Stouffer's method for combining z scores
        cv_zscore_est = sum(cv_zscore_folds) / sqrt(params$n_folds)
        cv_zscore_pval = 2*pnorm(abs(cv_zscore_est), lower.tail = FALSE)
        cv_pval_est = pchisq(-2 * sum(log(cv_pval_folds)), 2*params$n_folds, lower.tail = FALSE)
        
        if(fit$nzero[best_lam_ind] > 0) {
      
            # Number of predictors
            n_fit = fit$nzero[best_lam_ind]
            
            # Correlation between predicted and observed expression values for all samples ( as EpiXcan)
            y_all = predict(fit, as.matrix(train_data$predictors), s = 'lambda.min')
           
            corr_R = ifelse(sd(y_all) != 0, cor(y_all, train_data$response), 0)
            pval=cor.test(y_all, adj_expression)$p.value
            corr_R2 = corr_R**2
            
            # Adjusted correlation R2
            adj_R2 = 1 - (1 - corr_R2) * (nrow(train_data$predictors) - 1) / (nrow(train_data$predictors) - 1 - n_fit)
            
            # Get model coefficients for weights
            weights = fit$glmnet.fit$beta[which(fit$glmnet.fit$beta[,best_lam_ind] != 0), best_lam_ind]
            weighted_genes = names(fit$glmnet.fit$beta[,best_lam_ind])[which(fit$glmnet.fit$beta[,best_lam_ind] != 0)]
            
            # Update predictor weights
            train_data$predictor_weights = rep(0, ncol(train_data$predictors))
            names(train_data$predictor_weights) = colnames(train_data$predictors)
            train_data$predictor_weights[weighted_genes] = weights
            
            # Calculate additional metrics
            rmse = sqrt(mean((train_data$response - y_all)^2))
            mae = mean(abs(train_data$response - y_all))
           
            
        } else {
            # No predictors selected
            adj_R2 = NA
            corr_R = NA
            corr_R2 = NA
            rmse = NA
            mae = NA
            pval = NA
            train_data$predictor_weights = rep(0, ncol(train_data$predictors))
            names(train_data$predictor_weights) = colnames(train_data$predictors)
        }
        
    } else {
        # Model failed
        cv_R2_avg = NA
        cv_R2_sd = NA
        training_R2 = NA
        cv_rho_avg = NA
        cv_rho_se = NA
        cv_rho_avg_squared = NA
        cv_zscore_est = NA
        cv_zscore_pval = NA
        cv_pval_est = NA
        adj_R2 = NA
        corr_R = NA
        corr_R2 = NA
        rmse = NA
        mae = NA
        pval = NA
        train_data$predictor_weights = rep(0, ncol(train_data$predictors))
        names(train_data$predictor_weights) = colnames(train_data$predictors)
    }
    
    # Combine all metrics
    return(c(nested_perf, list(predictor_weights = train_data$predictor_weights,
        training_R2 = training_R2,
        cor_gene_all_data_pred = corr_R,
        adj_rsq_gene = adj_R2,
        rmse_avg = rmse,
        mae_avg = mae,
        pval_gene= pval,
        cv_R2_gene_avg = cv_R2_avg,
        cv_R2_gene_sd = cv_R2_sd,
        cv_rho_gene_avg = cv_rho_avg,
        cv_rho_gene_se = cv_rho_se,
        cv_rho_gene_avg_squared = cv_rho_avg_squared,
        cv_zscore_est_gene = cv_zscore_est,
        cv_zscore_pval_gene = cv_zscore_pval,
        cv_pval_est_gene = cv_pval_est,
        lambda_min = if(length(fit) > 0) fit$lambda.min else NA
    )))
}

#' Create comprehensive model summary matching CIS/EpiXcan format
#' @param gene Gene ID
#' @param train_data Training data list
#' @param perf Performance metrics
#' @param net_name Network name
#' @param mod_name Module name
#' @param params Model parameters
#' @return Data frame with comprehensive model summary
create_comprehensive_model_summary = function(gene, train_data, perf, net_name, mod_name, params) {
    gene_info = get_gene_info(train_data$annotations, gene)
    
    data.frame(
        gene_id = gene,
        gene_name = gene_info$name,
        gene_type = gene_info$type,
        chromosome = gene_info$chr,
        alpha = params$alpha,
        n_predictor_genes = ncol(train_data$predictors),  # Number of predictor genes
        n_genes_in_model = sum(train_data$predictor_weights != 0),  # Non-zero weights
        best_lambda = perf$lambda_min,
        inner_R2_avg = perf$R2_avg,
        inner_R2_sd = perf$R2_sd,
        inner_pval_est = perf$pval_est,
        inner_rho_avg = perf$rho_avg,
        inner_rho_se = perf$rho_se,
        inner_rho_zscore = perf$rho_zscore,
        inner_rho_avg_squared = perf$rho_avg_squared,
        inner_zscore_pval = perf$zscore_pval,
        cv_R2_gene_avg = perf$R2_avg,
        cv_R2_gene_sd = perf$R2_sd,
        cv_rho_gene_avg = perf$rho_avg,
        training_R2 = perf$training_R2,
        cv_rho_gene_se = perf$rho_se,
        cv_rho_gene_avg_squared = perf$rho_avg_squared,
        cv_zscore_est_gene = perf$rho_zscore,
        cv_zscore_pval_gene = perf$zscore_pval,
        cv_pval_est_gene = perf$pval_est,
        cor_gene_all_data_pred = perf$cor_gene_all_data_pred,
        adj_rsq_gene_all_data = perf$adj_rsq_gene_all_data,
        rmse_avg = perf$rmse_avg,
        pval_gene_all_data = perf$pval_gene_all_data,
        cv_adj_r2_gene_avg = perf$adj_rsq_gene_all_data,
        cv_adj_r2_gene_sd = perf$R2_sd,
        cv_pval_lm_gene_avg = perf$pval_gene_all_data,
        network = net_name,
        module = mod_name,
        stringsAsFactors = FALSE
    )
}

#' Create empty model summary for failed models
#' @param gene Gene ID
#' @param net_name Network name
#' @param mod_name Module name
#' @return Data frame with empty model summary
create_empty_model_summary = function(gene, net_name, mod_name) {
    data.frame(
        gene_id = gene,
        gene_name = NA,
        gene_type = NA,
        chromosome = NA,
        alpha = 0.5,
        n_predictor_genes = 0,
        n_genes_in_model = 0,
        best_lambda = NA,
        inner_R2_avg = NA,
        inner_R2_sd = NA,
        inner_pval_est = NA,
        inner_rho_avg = NA,
        inner_rho_se = NA,
        inner_rho_zscore = NA,
        inner_rho_avg_squared = NA,
        inner_zscore_pval = NA,
        cv_R2_gene_avg = NA,
        cv_R2_gene_sd = NA,
        cv_rho_gene_avg = NA,
        training_R2 = NA,
        cv_rho_gene_se = NA,
        cv_rho_gene_avg_squared = NA,
        cv_zscore_est_gene = NA,
        cv_zscore_pval_gene = NA,
        cv_pval_est_gene = NA,
        cor_gene_all_data_pred = NA,
        adj_rsq_gene_all_data = NA,
        rmse_avg = NA,
        pval_gene_all_data = NA,
        cv_adj_r2_gene_avg = NA,
        cv_adj_r2_gene_sd = NA,
        cv_pval_lm_gene_avg = NA,
        network = net_name,
        module = mod_name,
        stringsAsFactors = FALSE
    )
}

#' Extract comprehensive model weights matching CIS/EpiXcan format
#' @param predictor_weights Vector of predictor weights
#' @param train_data Training data list
#' @param gene Target gene
#' @param net_name Network name
#' @param mod_name Module name
#' @param gene_info Gene information
#' @return Data frame with comprehensive model weights
extract_comprehensive_model_weights = function(predictor_weights, train_data, gene, net_name, mod_name, gene_info) {
    if(ncol(train_data$predictors) == 0) return(data.frame())
    
    # Get predictor gene information
    predictor_genes = colnames(train_data$predictors)
    
    # Create weights data frame - one row per predictor gene
    weights_df = data.frame(
        gene           = rep(gene, length(predictor_genes)),
        gene_name      = rep(gene_info$name, length(predictor_genes)),
        gene_start     = rep(gene_info$start, length(predictor_genes)),
        gene_end       = rep(gene_info$end, length(predictor_genes)),
        strand         = rep(gene_info$strand, length(predictor_genes)),
        chr            = rep(gene_info$chr, length(predictor_genes)),
        predictor_gene = predictor_genes,  # One predictor gene per row
        weight         = predictor_weights,  # Use the passed-in weights
        alpha          = rep(0.5, length(predictor_genes)),
        network        = rep(net_name, length(predictor_genes)),
        module         = rep(mod_name, length(predictor_genes)),
        stringsAsFactors = FALSE
    )
    
    # Remove zero weights
    weights_df = weights_df[weights_df$weight != 0, ]
    
    return(weights_df)
}

#' Save consolidated results for a network
#' @param summaries Model summaries data frame
#' @param weights Model weights data frame
#' @param output_files List of output file paths
save_network_results = function(summaries, weights, output_files) {
    # Save summaries with error handling
    tryCatch({
        if(nrow(summaries) > 0) {
            write.table(summaries, file=output_files$summary, 
                        sep='\t', quote=FALSE, row.names=FALSE)
            message(sprintf("Saved %d model summaries to %s", 
                           nrow(summaries), output_files$summary))
        } else {
            message("No summaries to save")
        }
    }, error = function(e) {
        message(sprintf("Error saving summaries: %s", e$message))
    })
    
    # Save weights with error handling
    tryCatch({
        if(nrow(weights) > 0) {
            write.table(weights, file=output_files$weights,
                        sep='\t', quote=FALSE, row.names=FALSE)
            message(sprintf("Saved %d model weights to %s", 
                           nrow(weights), output_files$weights))
        } else {
            message("No weights to save")
        }
    }, error = function(e) {
        message(sprintf("Error saving weights: %s", e$message))
    })
}

#' Filter module list to remove unwanted modules (e.g. grey module)
#' @param moduleList List of network modules
#' @return Filtered module list
filter_module_list = function(moduleList) {
    moduleList = moduleList[!names(moduleList) %in% c("Fromer2016_case","Pergola2017","hippo.noQSVAremoved","dentate.QSVAremoved")]
    map(names(moduleList) %>% set_names(.,.), function(net) {
        moduleList[[net]][!grepl("grey", names(moduleList[[net]]))]
    })
}

#' Prepare training data for a single gene
#' @param gene Gene ID
#' @param module_data Module data list
#' @return List of prepared training data
prepare_gene_data = function(gene, module_data) {
    tryCatch({
        # Remove target gene from predictors if present
        if(gene %in% colnames(module_data$imputed)) {
            predictors = module_data$imputed[,-which(colnames(module_data$imputed) == gene)]
        } else {
            predictors = module_data$imputed
        }
        
        # Get response variable (target gene expression)
        if(gene %in% colnames(module_data$real)) {
            response = module_data$real[,which(colnames(module_data$real) == gene),drop=F]
        } else {
            return(NULL)
        }
        
        # Match samples between predictors and response
        common_samples = intersect(rownames(predictors), rownames(response))
        if(length(common_samples) == 0) {
            message(sprintf("No common samples found for gene %s", gene))
            return(NULL)
        }
        
        predictors = predictors[common_samples,,drop=F]
        response = response[common_samples,,drop=F]
        
        # Scale response variable
        response = scale(response, center = TRUE, scale = TRUE)
        
        # Get predictor weights (will be calculated during model training)
        predictor_weights = rep(0, ncol(predictors))
        names(predictor_weights) = colnames(predictors)
        
        list(
            predictors = predictors,
            response = response,
            predictor_weights = predictor_weights,
            annotations = module_data$annotations
        )
        
    }, error = function(e) {
        message(sprintf("Error preparing data for gene %s: %s", gene, e$message))
        return(NULL)
    })
}

#' Get file path for different file types
#' @param file_type Type of file
#' @param brain_region Brain region name
#' @param net_name Network name
#' @param base_dir Base directory
#' @return File path
get_file_path = function(file_type, brain_region, net_name, base_dir) {
    patterns = list(
        summary = file.path(base_dir, "summary", brain_region, 
                           paste0(net_name, "_summary.txt")),
        weights = file.path(base_dir, "weights", brain_region, 
                           paste0(net_name, "_weights.txt")),
        predictions = file.path(base_dir, "predictions", brain_region, 
                               paste0(net_name, "_predictions.txt"))
    )
    
    patterns[[file_type]]
} 