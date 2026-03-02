#####################################################################
# INGENE Prediction Helper Functions
#####################################################################
# Helper functions for applying trained INGENE models to new data and 
# averaging predictions across networks
#####################################################################

#' Compute predictions from trained weights for a single network
#' @param test.pred Test data frame with expression data
#' @param mod.weights Model weights from training
#' @param network Network name
#' @param output_file Full path to output file
#' @return Prediction matrix or NA if error/already exists
computePredFromWeights <- function(test.pred, mod.weights, network, output_file, n_cores) {
    tryCatch({
        # Input validation
        if (is.null(test.pred) || nrow(test.pred) == 0) {
            warning("Empty or NULL test data provided")
            return(NULL)
        }
        if (is.null(mod.weights) || nrow(mod.weights) == 0) {
            warning("Empty or NULL model weights provided")
            return(NULL)
        }
        if (!file.exists(dirname(output_file))) {
            dir.create(dirname(output_file), recursive = TRUE)
        }
        
        # Remove 0-weight entries to optimize computation
        mod.weights = mod.weights[mod.weights$weight != 0, ]
        if(nrow(mod.weights) == 0){
            warning("No non-zero weights found in model")
            return(NULL)
        }
        
        # Filter predictors present in test data and group by target genes
        mod.weights_sub = mod.weights %>%
            group_by(gene) %>%
            filter(predictor_gene %in% colnames(test.pred))
        
        # Get list of genes we can predict
        predictable_genes = unique(mod.weights_sub$gene)
        if(length(predictable_genes) == 0) {
            warning("No predictable genes found in the data")
            return(NULL)
        }
        
        # Filter weights to specific network
        mod.weights_sub = mod.weights_sub[which(mod.weights_sub$network == network),]
        if(nrow(mod.weights_sub) == 0) {
            warning(paste("No weights found for network:", network))
            return(NULL)
        }
        
        # Check if predictions already exist
        if(file.exists(output_file)) {
           return( message("File exists"))
            
        }
        
        # Get common predictors between test data and model
        samples = test.pred$IID
        common_predictors = intersect(colnames(test.pred), mod.weights_sub$predictor_gene)
        if(length(common_predictors) == 0) {
            warning("No common predictors found between test data and model")
            return(NULL)
        }
        
        # Extract relevant predictors from test data
        test.predictors = test.pred %>%
          dplyr::select(IID, all_of(common_predictors))
        
        # Compute predictions in parallel
        module_predictions_list = mclapply(predictable_genes, function(gene) {
          message(sprintf("Gene %d/%d", which(predictable_genes==gene), length(predictable_genes)))
          
          test.partners = test.predictors %>%
            dplyr::select(matches(mod.weights_sub$predictor_gene[
              mod.weights_sub$gene == gene])) %>%
            as.data.frame()
          
          weights = mod.weights_sub[mod.weights_sub$gene == gene, ] %>%
            filter(predictor_gene %in% colnames(test.partners)) %>%
            dplyr::select(gene, weight, predictor_gene) %>%
            dplyr::rename(target = "gene", predictor = "predictor_gene")
          
          weights$weight = as.numeric(weights$weight)
          
          if (!identical(colnames(test.partners), weights$predictor)) {
            warning(paste("Predictor mismatch for gene:", gene))
            return(tibble(!!gene := rep(NA, length(samples))))
          }
          
         
          res.mod = matrix(NA, nrow = length(samples), ncol = length(weights$predictor))
          for (col in seq_along(weights$predictor)) {
            res.mod[, col] = test.partners[, col] * weights$weight[col]
          }
          
          res.gene = rowSums(res.mod)
          
          return(tibble(!!gene := res.gene))
        }, mc.cores = n_cores, mc.preschedule = FALSE)
        
        # Combine predictions
        module_predictions = bind_cols(module_predictions_list)
        
        # Combine with sample IDs
        module_res = bind_cols(tibble(IID = samples), module_predictions)
        
        # Validate final output
        if(ncol(module_res) <= 1) {
            warning("No predictions generated")
            return(NULL)
        }
        
        # Save predictions
        tryCatch({
            write.table(module_res, file = output_file, 
                       sep = "\t", quote = FALSE, row.names = FALSE)
        }, error = function(e) {
            warning("Failed to save predictions to file: ", e$message)
        })
        
        return(module_res)
        
    }, error = function(e) {
        warning("Error computing predictions for network ", network, ": ", e$message)
        return(NULL)
    })
}


