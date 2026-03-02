#####################################################################
# Helper Functions for cis-eQTL Analysis
#####################################################################

#' Generate stratified fold IDs based on diagnosis
#' @description Creates cross-validation fold assignments while maintaining balanced 
#'             representation of diagnosis groups across folds
#' @param n_folds Number of folds for cross-validation
#' @param covariates Data frame containing sample covariates including diagnosis (Dx)
#' @return Vector of fold assignments (1 to n_folds) for each sample
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

#' Filter genes and SNPs within cis window
#' @description Identifies SNP-gene pairs within a specified genomic distance (cis window)
#' @details For each gene, finds all SNPs that lie within the specified window size
#'          upstream or downstream of the gene boundaries
#' @param gene_annot Data frame with gene annotations (chr, tss positions)
#' @param snp_annot Data frame with SNP annotations (chromosome, position)
#' @param chromosome Current chromosome being analyzed
#' @param cis_window Size of cis window in base pairs
#' @return List containing filtered genes, SNPs and their pairwise associations
filter_cis_pairs = function(gene_annot, snp_annot, chromosome, cis_window) {
    # Filter genes and SNPs on current chromosome
    chr_genes = gene_annot[gene_annot$chr == chromosome, ]
    chr_snps  = snp_annot[snp_annot$Chromosome == chromosome, ]
    
    if(nrow(chr_genes) == 0 || nrow(chr_snps) == 0) return(NULL)
    
    # Create pairs within cis window
    cis_pairs = list()
    for(i in 1:nrow(chr_genes)) {
        gene_id = chr_genes$gencodeID[i]
        gene_tss = chr_genes[i,"tss"]
        gene_end = chr_genes[i,"end"] + cis_window
        
        # Find SNPs in cis window for this gene
        cis_snps = chr_snps[chr_snps$Position >= gene_tss - cis_window & 
                           chr_snps$Position <= gene_tss + cis_window,]
        
        if(nrow(cis_snps) > 0) {
            cis_pairs[[gene_id]] = unique(cis_snps$SNP)
        }
    }
    
    return(list(
        genes = chr_genes,
        snps = chr_snps,
        pairs = cis_pairs
    ))
}

# Function for fixed-effect inverse variance meta-analysis
meta_effect = function(betas, ses) {
  weights = 1 / (ses^2)
  beta_meta = sum(weights * betas) / sum(weights)
  se_meta = sqrt(1 / sum(weights))
  z_meta = beta_meta / se_meta
  p_meta = 2 * pnorm(-abs(z_meta))
  return(c(beta_meta, se_meta, p_meta))
}

#' Run robust linear model analysis for cis-eQTLs with cross-validation
#' @description Performs cross-validated robust linear regression to identify eQTLs
#' @details For each SNP-gene pair:
#'          1. Splits data into train/test folds
#'          2. Fits robust linear models controlling for covariates
#'          3. Calculates association statistics and rankings
#' @param gene Gene identifier being analyzed
#' @param expression Expression data matrix for the gene
#' @param genotypes Matrix of genotypes for cis-window SNPs
#' @param snp_annot SNP annotation information
#' @param covariates Sample covariates to include in the model
#' @param n_folds Number of cross-validation folds
#' @param n_cores Number of CPU cores for parallel processing
#' @param output_dir Directory to save results
#' @param data_dir Directory containing input data
#' @param cv_fold_id Identifier for CV fold assignments
#' @param brain_region Brain region being analyzed
#' @return Data frame of significant eQTL associations
run_rlm_analysis = function(gene, expression, genotypes, snp_annot, covariates, n_folds=4, n_cores = 1,
                           output_dir, data_dir, cv_fold_id, brain.region) {
  
    
    # Match samples between expression, genotype and covariates
    expression = expression[which(rownames(expression) %in% rownames(genotypes)),,drop=F]
    genotypes  = genotypes[which(rownames(genotypes) %in% rownames(expression)),,drop=F]
    expression = expression[match(rownames(genotypes), rownames(expression)),,drop=F]
    stopifnot(identical(rownames(expression), rownames(genotypes)))
    
    covariates = covariates[which(rownames(covariates) %in% rownames(genotypes)),,drop=F]
    covariates = covariates[match(rownames(genotypes), rownames(covariates)),,drop=F]
    stopifnot(identical(rownames(covariates), rownames(genotypes)))
    
   
    ## Remove SNPs with missing values (Elastic Net training does not handle NAs)
    genotypes = genotypes[, apply(genotypes, 2, function(x) all(!is.na(x)))]
    
    
    # Generate fold assignments
    fold_file = file.path(data_dir, 
                          paste0(cv_fold_id, '_cv_4fold_', brain.region, '_ids.RData'))
    
    fold_ids = if(!file.exists(fold_file)) {
      ids    = generate_fold_ids_stratified(n_folds, covariates)
      save(ids, file = fold_file)
      ids
    } else {
      get(load(fold_file))
    }
    
    
    # Run cross-validation
    cv_results = lapply(1:n_folds, function(fold) {
        # Split into train/test sets
        train_idx = ids != fold
        test_idx  = ids == fold
        
        x_train    = genotypes[train_idx, ]
        y_train    = expression[train_idx, ]
        covs_train = covariates[train_idx, ]
        
        ## Run rlm for each SNP
        QTL.rlm = mclapply(colnames(x_train), function(snp) {
          if (length(unique(x_train[, snp])) == 1) return(NULL)  # Skip SNPs with no variability
          
          df = data.frame(y_train = y_train, snp_value = x_train[, snp])  # Base dataframe
          
          # Include covariates if available
          #covar_terms = colnames(covs_train)[grepl("Age|Sex|Dx|C[1-5]|PC", colnames(covs_train))]
          #df          = cbind(df, covs_train)  # Add covariates to the dataframe
          # df$Sex      = as.numeric(df$Sex)
          # df$Dx       = as.numeric(df$Dx)
          #formula_str = paste("y_train ~", paste(c(covar_terms, "snp_value"), collapse = " + "))
          formula_str = paste("y_train ~ snp_value")
        
          form = as.formula(formula_str)
          
          # Fit RLM model with error handling
          tryCatch(rlm(form, data = df), error = function(e) NULL)
        }, mc.cores = n_cores, mc.preschedule = FALSE)
        
        
        names(QTL.rlm) = colnames(x_train)
        QTL.rlm        = QTL.rlm[!sapply(QTL.rlm, is.null)]
        
        # Extract coefficients and p-values
        coef_stats = lapply(QTL.rlm, function(model) {
          coef_summary = summary(model)$coefficients  # Extract coefficients
          snp_idx       = which(rownames(coef_summary) == "snp_value")  # Find SNP row
          coef_summary[snp_idx, ]
        })
        
        p_values  = lapply(QTL.rlm, function(model) {
          f.robftest(model, var = "snp_value")$p.value
        })
        
        # Format results for this fold
        fold_results = data.frame(
          snp     = names(QTL.rlm),
          beta    = unlist(lapply(coef_stats, function(x) x[1])),
          st_err  = unlist(lapply(coef_stats, function(x) x[2])),
          t_value = unlist(lapply(coef_stats, function(x) x[3])),
          p_value = unlist(p_values)
        )
        
        # Add fold-specific rankings
        fold_results = fold_results[order(fold_results$p_value),]
        fold_results$rank = rank(fold_results$p_value)
        
        return(fold_results)
    })
    
    # Get significant SNPs and save results
    significant_eqtls = getSignificantQTLs(
        cv_list     = cv_results,
        snp_pos     = snp_annot, 
        genotype    = genotypes
    )
    
    significant_eqtls$gene = gene
    
    return(significant_eqtls)
}


## Function to perform pruning
#' Calculate linkage disequilibrium between SNPs
#' @description Computes r-squared values between a target SNP and neighboring SNPs
#'             within a specified genomic window
#' @details For a given SNP:
#'          1. Identifies all SNPs within the flanking window
#'          2. Calculates pairwise r-squared values
#'          3. Returns SNPs exceeding the LD threshold for pruning
#' @param genotypes Matrix of genotype values
#' @param flank Size of flanking region (in base pairs) to check for LD
#' @param association Data frame containing SNP positions and association statistics
#' @param snp Index of the target SNP in the association data frame
#' @param rsquared R-squared threshold for considering SNPs in LD
#' @return List containing:
#'         - cut: r-squared values above threshold
#'         - new.genotypes: Genotype matrix with LD SNPs removed
#'         - new.association: Association results with LD SNPs removed
#'         - snp: IDs of SNPs removed due to LD
SNPr2 = function(genotypes, flank, association, snp, rsquared) {
  #print(snp)
  index = which(association[,2]==association[snp,2] &
                  association[,3]>= association[snp,3]-flank &
                  association[,3]<= association[snp,3]+flank)
  
  marker         = row.names(association)[index]
  window         = genotypes[,which(colnames(genotypes) %in% marker), drop=F]
  r2             = cor(window, genotypes[,association[snp,1]], use="pairwise.complete.obs")^2
  cut            = subset(r2, r2[,1] >= rsquared)
  snp_r2         = row.names(cut)[!(row.names(cut) %in% association[snp,1])]
  new_genotypes  = genotypes[,!(colnames(genotypes) %in% snp_r2)]
  new_assocation = association[!(association[,1] %in% snp_r2),]
  results        = list(cut=cut,
                    new.genotypes=new_genotypes,
                    new.association=new_assocation,
                    snp=snp_r2)
  return(results)
}

#' Prune SNPs in linkage disequilibrium
#' @description Removes highly correlated SNPs based on LD structure to identify
#'             independent association signals
#' @details Iterative process:
#'          1. Takes SNPs in order of association strength
#'          2. For each SNP, identifies and removes others in high LD
#'          3. Updates association statistics after pruning
#'          4. Applies multiple testing corrections
#' @param annot Data frame with SNP annotations (SNP ID, chromosome, position)
#' @param eqtls Data frame containing eQTL association results
#' @param x Matrix of genotype values for all SNPs
#' @return Data frame of pruned eQTL associations with updated statistics:
#'         - Original association metrics
#'         - FDR and Bonferroni corrected p-values
#'         - Updated rankings after pruning
ciseQTL_pruning = function(annot, eqtls, x) {
  
  # Filter and match annot based on eqtls
  annot = annot[annot$SNP %in% rownames(eqtls), ]
  annot = annot[match(rownames(eqtls), annot$SNP), ]
  
  stopifnot(identical(rownames(eqtls), annot$SNP))
  
  # Create association data frame
  assoc = cbind(annot[, c("SNP", "Chromosome", "Position")], eqtls)
  colnames(assoc)[1] ="Marker"
  assoc              = assoc[order(assoc$meta_p), ]
  rownames(assoc)    = assoc$Marker
  
  # Pruning parameters
  f = 250000
  rsq = 0.9
  
  # Use a while loop for pruning
  i = 1
  while (i <= nrow(assoc)) {
    pruning = SNPr2(genotypes = x, flank = f, association = assoc, snp = i, rsquared = rsq)
    x       = pruning$new.genotypes
    assoc   = pruning$new.association
    i        = i + 1
  }
  
  # Adjust p-values and ranks
  assoc$FDR        = p.adjust(assoc$new_pval, method = "fdr")
  assoc$Bonferroni = p.adjust(assoc$new_pval, method = "bonferroni")
  assoc$rank       = rank(assoc$new_pval)
  
  return(assoc)
}


#' Get significant SNPs across CV folds
#' @description Processes cross-validation results to identify consistent eQTLs
#' @details 1. Finds SNPs present in all folds
#'          2. Combines statistics across folds
#'          3. Performs LD pruning on final results
#' @param gene Gene identifier
#' @param cv_list List of results from each CV fold
#' @param snp_pos SNP position information
#' @param genotype Full genotype data
#' @return Data frame of significant pruned eQTL associations
getSignificantQTLs = function(gene, cv_list, snp_pos, genotype) {
     
      
      # Get common SNPs across folds
      common_snps = Reduce(intersect, lapply(cv_list, function(x) x$snp))
      
      # Filter each fold to retain only common SNPs
      cv_list = lapply(cv_list, function(fold) {
        fold  = fold[fold$snp %in% common_snps, , drop = FALSE]
        return(fold)
      })
      
      # Ensure SNP order is the same across all folds (based on the first fold)
      snp_order = cv_list[[1]]$snp
      
      # Reorder each fold to match the first fold SNP order
      cv_list   = lapply(cv_list, function(fold) {
        fold = fold[match(snp_order, fold$snp), , drop = FALSE]
        return(fold)
      })
      
      # Merge results across folds
      merged_df = Reduce(cbind, lapply(cv_list, function(x) x))
      
      # Extract relevant statistics
      rank_cols = merged_df[, grepl('rank', colnames(merged_df)), drop = FALSE]
      pval_cols = merged_df[, grepl('p_value', colnames(merged_df)), drop = FALSE]
      stat_cols = merged_df[, grepl('snp|beta|st_err|t_value', colnames(merged_df)), drop = FALSE]
      
      # Compute combined statistics
      merged_results = data.frame(
        new_rank = rank(rowProds(as.matrix(rank_cols))),
        new_pval = rowProds(as.matrix(pval_cols))
      )
      
      # Extract beta and se columns for meta-analysis
      beta_cols = merged_df[, grepl('^beta(\\.|$)', colnames(merged_df)), drop = FALSE]
      se_cols   = merged_df[, grepl('^st_err(\\.|$)', colnames(merged_df)), drop = FALSE]
      
      # Apply meta-analysis row-wise
      meta_results = t(sapply(seq_len(nrow(beta_cols)), function(i) {
        betas = as.numeric(beta_cols[i, ])
        ses   = as.numeric(se_cols[i, ])
        meta_effect(betas, ses)
      }))
      colnames(meta_results) = c("meta_beta", "meta_se", "meta_p")
      
      # Combine meta-analysis results with merged_df
      final_results = cbind(merged_results, meta_results, rank_cols, pval_cols, stat_cols)
      
      # Sort by meta_p
      final_results = final_results[order(final_results$meta_p), ]
      rownames(final_results) = final_results$snp
      
      # Prune SNPs in LD
      pruned_results = ciseQTL_pruning(annot = snp_pos, eqtls = final_results, x = genotype)
      
      return(pruned_results)
    }
    
   


