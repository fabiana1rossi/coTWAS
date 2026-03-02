#####################################################################
# INGENE Network Prediction Averaging Helper Functions
#####################################################################
# This script contains functions for averaging and validating predictions
# across different networks in the INGENE pipeline. The main function
# 'ensembleAcrossNet' performs three key steps:
# 1. Validates predictions against true expression data
# 2. Filters genes based on performance metrics
# 3. Averages predictions for genes present in multiple networks
# 4. Preserves unique gene predictions from single networks
#####################################################################

#####################
## Average predictions across networks
########################

#' Ensemble predictions across multiple networks with validation
#' @param brain.region Brain region being analyzed
#' @param model Model type being used (cis/epixcan for INGENE, "MODULE" for MODULE)
#' @param net.names Names of networks to ensemble
#' @param true.expr True expression data for validation
#' @param outputdir Directory to save results
#' @param outname Output file name suffix
#' @param network_predictions List of prediction dataframes from different networks
#' @param pipeline Name of the pipeline ("INGENE" or "MODULE")
#' @return Dataframe with averaged predictions across networks
ensembleAcrossNet = function(brain.region, model, net.names, true.expr, outputdir, outname, testing_geno, network_predictions, pipeline, genestokeep=NULL, n_cores) {
  # Load required base packages
  library(parallel)
  
  ##########################
  ### STEP 1: Initialize and Clean Network Predictions
  ##########################
  pred.expr_lst = lapply(network_predictions, function(x) {
    if (nrow(x) != 0) return(x)
  })
  names(pred.expr_lst) = net.names
  pred.expr_lst = pred.expr_lst[!sapply(pred.expr_lst, is.null)]
  net.names = names(pred.expr_lst)
  
  ############
  ### STEP 2: Validate Predictions Against True Expression
  ############
  pred.expr_lst = mclapply(net.names, function(net) {
    cat(net, "\n")
    
    if (!is.null(genestokeep)) {
      keep_genes = c("IID", "FID", genestokeep$gene_name[genestokeep$network == net])
      pred.expr_lst[[net]] = pred.expr_lst[[net]][, colnames(pred.expr_lst[[net]]) %in% keep_genes, drop = FALSE]
    }
    
    testing_perf = run_model_analysis(real.expr = true.expr,
                                      predicted.expr = pred.expr_lst[[net]],
                                      cov_df = NULL,
                                      threshold = 'yes')
    
    tmp = pred.expr_lst[[net]]
    if (!is.null(testing_perf) && nrow(testing_perf) > 0) {
      tmp = tmp[, colnames(tmp) %in% c('IID', 'FID', testing_perf$gene), drop = FALSE]
      return(list(testing_perf, tmp))
    } else {
      return(NA)
    }
  }, mc.cores = n_cores)
  
  names(pred.expr_lst) = net.names
  pred.expr_lst = pred.expr_lst[!sapply(pred.expr_lst, function(x) all(is.na(x)))]
  
  if (length(pred.expr_lst) == 0) {
    message("No significant genes across networks")
    return(NULL)
  }
  
  # Combine performance metrics
  testing_perf = do.call(rbind, lapply(names(pred.expr_lst), function(net) {
    tmp = pred.expr_lst[[net]][[1]]
    tmp$network = net
    tmp$region = brain.region
    return(tmp)
  }))
  
  # Save performance
  if (!is.null(testing_perf)) {
    perf_file = if (pipeline == "INGENE") {
      paste0(outputdir, "/INGENE_", testing_geno, "_performance_", brain.region, ".RData")
    } else {
      paste0(outputdir, "/MODULE_", testing_geno, "_performance_", brain.region, ".RData")
    }
    save(testing_perf, file = perf_file)
  }
  
  pred.expr_lst = setNames(lapply(names(pred.expr_lst), function(net) {
    pred.expr_lst[[net]][[2]]
  }), names(pred.expr_lst))
  
  all_iids = lapply(pred.expr_lst, function(x) x$IID)
  
  # Check for valid IIDs in each network
  valid_networks = sapply(all_iids, function(iids) {
    # Check if IIDs exist and are not all NA/empty
    if (is.null(iids) || length(iids) == 0) {
      return(FALSE)
    }
    # Check if all IIDs are valid (not NA, not empty string)
    all(!is.na(iids) & iids != "" & !is.null(iids))
  })
  
  # Remove networks with invalid IIDs
  if (!all(valid_networks)) {
    invalid_networks = names(valid_networks)[!valid_networks]
    message("Warning: Discarding ", length(invalid_networks), " network(s) with invalid IIDs: ", 
            paste(invalid_networks, collapse = ", "))
    
    pred.expr_lst = pred.expr_lst[valid_networks]
    net.names = names(pred.expr_lst)
    
    if (length(pred.expr_lst) == 0) {
      stop("No networks with valid IIDs remaining after filtering")
    }
    
    message("Remaining networks with valid IIDs: ", paste(net.names, collapse = ", "))
  }
  
  # Now check if all remaining networks have the same IIDs
  all_iids = lapply(pred.expr_lst, function(x) x$IID)
  
  if (!all(sapply(all_iids, function(x) identical(x, all_iids[[1]])))) {
    stop("All prediction dataframes must have the same IIDs. Networks with different IIDs will be discarded.")
  }
  
  sample_ids = all_iids[[1]]
  message("All networks have ", length(sample_ids), " samples with matching IIDs")
  
  ##########################
  ### STEP 3: Average Predictions Across Networks
  ##########################
  
  genes_by_network = lapply(pred.expr_lst, colnames)
  
  gene_list = unlist(genes_by_network, use.names = FALSE)
  network_list = rep(names(genes_by_network), times = sapply(genes_by_network, length))
  
  gene_table = data.frame(gene = gene_list, network = network_list, stringsAsFactors = FALSE)
  gene_table = gene_table[!gene_table$gene %in% c("IID", "FID"), , drop = FALSE]
  
  gene_summary = aggregate(network ~ gene, data = gene_table, 
                           FUN = function(x) list(unique(x)))
  gene_summary$n_networks = sapply(gene_summary$network, function(x) length(x))
  
  common_genes = gene_summary$gene[gene_summary$n_networks > 1]
  unique_genes = gene_summary$gene[gene_summary$n_networks == 1]
  
  # Common gene averaging
  common_genes_predictions = lapply(common_genes, function(gene) {
    networks = unlist(gene_summary$network[gene_summary$gene == gene])
    
    # Get the maximum number of rows across all networks for this gene
    max_rows = max(sapply(pred.expr_lst[networks], function(df) nrow(df)))
    
    # Filter networks to only include those with the maximum number of rows
    valid_networks = networks[sapply(pred.expr_lst[networks], function(df) nrow(df) == max_rows)]
    
    # Log if networks were excluded due to fewer rows
    if (length(valid_networks) < length(networks)) {
      excluded_networks = setdiff(networks, valid_networks)
      message("Gene ", gene, ": Excluded ", length(excluded_networks), " network(s) with fewer rows: ", 
              paste(excluded_networks, collapse = ", "))
    }
    
    # If only one network remains after filtering, use it directly
    if (length(valid_networks) == 1) {
      out = data.frame(pred.expr_lst[[valid_networks[1]]][, gene, drop = FALSE])
      rownames(out) = sample_ids
      colnames(out) = gene
      return(out)
    }
    
    # Otherwise, average across the valid networks
    gene_preds = do.call(cbind, lapply(pred.expr_lst[valid_networks], function(df) df[, gene, drop = FALSE]))
    avg_pred = rowMeans(gene_preds, na.rm = TRUE)
    out = data.frame(avg_pred)
    rownames(out) = sample_ids
    colnames(out) = gene
    return(out)
  })
  names(common_genes_predictions) = common_genes
  
  # Unique gene selection
  unique_genes_predictions = lapply(unique_genes, function(gene) {
    networks = unlist(gene_summary$network[gene_summary$gene == gene])
    pred.expr_lst[[networks[1]]][, gene, drop = FALSE]
  })
  names(unique_genes_predictions) = unique_genes
  
  # Final prediction matrix
  all_predictions = do.call(cbind, c(common_genes_predictions, unique_genes_predictions))
  all_predictions$IID = all_predictions$FID = sample_ids
  
  ensemble_perf = run_model_analysis(real.expr = true.expr,
                                     predicted.expr = all_predictions,
                                     cov_df = NULL,
                                     threshold = 'yes')
  
  if (!is.null(ensemble_perf) && nrow(ensemble_perf) > 0) {
    ensemble_perf$network = "Averaged_Ensemble"
    ensemble_perf$region = brain.region
    testing_perf = rbind(testing_perf, ensemble_perf)
    
    perf_file = if (pipeline == "INGENE") {
      paste0(outputdir, "/INGENE_", testing_geno, "_performance_", brain.region, ".RData")
    } else {
      paste0(outputdir, "/MODULE_", testing_geno, "_performance_", brain.region, ".RData")
    }
    save(testing_perf, file = perf_file)
    
    significant_genes = testing_perf$gene[testing_perf$network == "Averaged_Ensemble"]
    keep_cols = colnames(all_predictions) %in% c("IID", "FID", significant_genes)
    all_predictions = all_predictions[, keep_cols, drop = FALSE]
    
    output_file = if (pipeline == "INGENE") {
      paste0(outputdir, '/', brain.region, '_', model, '_AveragedNetworks_', outname, '.txt')
    } else {
      paste0(outputdir, '/', brain.region, '_AveragedNetworks_', outname, '.txt')
    }
    
    write.table(all_predictions,
                file = output_file,
                col.names = TRUE,
                row.names = TRUE,
                sep = '\t',
                quote = FALSE)
    return(all_predictions)
  } else {
    message("No significant genes across networks")
    return(NULL)
  }
}

#' Helper function to average predictions for a single gene
#' @param predictions Matrix of predictions to average
#' @return Vector of averaged predictions
#' @description Calculates the mean of predictions across networks,
#'              handling missing values appropriately
averagePredictions = function(predictions) {
    rowMeans(predictions, na.rm = TRUE)
}

