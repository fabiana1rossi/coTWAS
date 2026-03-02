#####################################################################
# Network Training Pipeline
#####################################################################
# This script processes the coeQTL results to train network models
# It performs the following steps:
# 1. Loads module/network-specific coeQTL results
# 2. Processes significant coeQTLs for each network
# 3. Trains and validates network models
# 4. Creates extra and weights files for each network
#####################################################################

#################################################
## Setup and Configuration
#################################################

# Load required libraries
required_packages = c(
    "dplyr", "purrr", "parallel", "glmnet",
    "magrittr", "moments",
    "tibble","limma","stringr"
)

# Install and load packages
for(pkg in required_packages) {
    if(!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
        library(pkg, character.only = TRUE)
    }
}


## set working directory
# setwd("")


# Source helper functions
source("./STEP3.Training_MODULE/training_helper_functions.R")



#################################################
## Configuration
#################################################
base_dir              = getwd()

CONFIG = list(
    data_file_path        = paste0(base_dir, "/data/"),
    rnaseq                = paste0(getwd(), "/rnaseq/your_rnaseq.RData"),
    network_dir           = paste0(getwd(), "/STEP3.Training_MODULE/output/"),
    output_dir            = paste0(getwd(), "/STEP3.Training_MODULE/output/"),
    brain_regions         = c('dlpfc'),
    n_cores               = 1
)


#################################################
## Main Analysis Functions
#################################################

#' Process a single brain region
#' @param brain_region Brain region name
#' @param config Configuration parameters
#' @return NULL (saves results to files)
process_brain_region = function(brain_region, config) {
    message("Processing brain region: ", brain_region)
    
    ## Create summary and weights folders to store training results
    summary_output_dir = file.path(config$output_dir, "summary",brain_region)
    if (!dir.exists(summary_output_dir)) {
        dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    weights_output_dir = file.path(config$output_dir, "weights",brain_region)
    if (!dir.exists(weights_output_dir)) {
        dir.create(weights_output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    # Load expression and covariate data
    training_data = get(load(config$rnaseq))
    expression    = as.data.frame(training_data[[brain_region]]$assays$expression)
    covariates    = as.data.frame(training_data[[brain_region]]$colData)
    gene_annot    = as.data.frame(training_data[[brain_region]]$rowData)
    
    
    # Get network files containing modules data (genotype, expression, coeQTLs)
    network_files  = list.files(
        config$network_dir,
        pattern    = paste0(brain_region, "_network"),
        full.names = TRUE
    )
    
    # Build a tibble with parsed network and module names
    file_tbl = tibble(file = network_files) %>%
      dplyr::mutate(
        filename = basename(file),
        network  = str_extract(filename, "^.*(?=_module_)"),
        module  = str_extract(filename, "(?<=_module_).*(?=\\.RData)")
      )
    
    # Nest into a named list of lists using purrr
    network_list = file_tbl %>%
      split(.$network) %>%
      map(~ set_names(.$file, .$module))
    
    
    # Process each network
    for(network_name in names(network_list)) {
        message("Processing network file: ", network_name)
        
      if(file.exists(file.path(weights_output_dir, paste0(network_name, ".weights.txt")))){message("Network ", network_name," already done");next}
      
       # Process each module in the network
        network_results = map(names(network_list[[network_name]]) %>% set_names(.,.), function(module_name) {
            # Get module specific file
            module_file  = network_files[grepl(sprintf("%s_module_%s.RData",network_name, module_name),network_files)]
            
            module = tryCatch({
              get(load(module_file))
            }, error = function(e) {
              message("Error loading module: ", e$message)
              NA
            })
            
            ## Get significant coeQTLs selected as the top 5% based on pval
            if(length(module) < 7){return(NA)} ## No eQTLs
            
            ## Get significant coeQTLs selected as the top 5% based on pval
            if(length(unique(module$coeqtls_pruned$Marker)) <=1){
              message("Not significant coeQTLs")
              return(NA)}
            
            module$coeqtls_pruned = module$coeqtls_pruned[order(module$coeqtls_pruned$rank,decreasing = FALSE), ]
            significant_coeqtls   = module$coeqtls_pruned$Marker[module$coeqtls_pruned$new_pval <= quantile(module$coeqtls_pruned$new_pval, 0.05)]
           
            ## Subset module genes that are present in the expression matrix and annotation
            module$genes  = module$genes[which(module$genes %in% intersect(colnames(expression), gene_annot$gencodeID))]
            
            # Train models for each gene in the module
            gene_models = mclapply(module$genes, function(gene) {
                print(paste0("Training model for gene: ", which(module$genes == gene), "/", length(module$genes)))
                
                tryCatch({
                    result = train_elastic_net_module(
                        brain.region        = brain_region,
                        network_name        = network_name,
                        module_name         = module_name,
                        data_file_path      = config$data_file_path,
                        significant_coeqtls = significant_coeqtls,
                        x                   = module$geno.mod,
                        y                   = expression[, which(colnames(expression) %in% module$genes), drop=FALSE],
                        gene                = gene,
                        gene_annot          = gene_annot,
                        snp_annot           = module$map.mod,
                        covs                = covariates,
                        cv_fold_id          = "MODULE"
                    )
                    
                    # Check if result contains valid model summary
                    if (!is.null(result$model_summary)) {
                        # Ensure all required statistics are present, if not set to NA
                        required_stats = c("cv_R2_pc1_avg", "cv_R2_pc1_sd", "cv_rho_pc1_avg", 
                                         "cv_rho_pc1_se", "cv_rho_pc1_avg_squared", "cv_zscore_est_pc1",
                                         "cv_zscore_pval_pc1", "cv_pval_est_pc1", "cv_R2_gene_avg",
                                         "cv_R2_gene_sd", "cv_rho_gene_avg", "cv_rho_gene_se",
                                         "cv_rho_gene_avg_squared", "cv_zscore_est_gene",
                                         "cv_zscore_pval_gene", "cv_pval_est_gene")
                        
                        for (stat in required_stats) {
                            if (is.null(result$model_summary[[stat]]) || 
                                is.na(result$model_summary[[stat]]) || 
                                !is.numeric(result$model_summary[[stat]])) {
                                result$model_summary[[stat]] = NA
                            }
                        }
                    }
                    
                    return(result)
                }, error = function(e) {
                    warning(paste("Error in model training for gene:", gene, "- Error:", e$message))
                    # Return a minimal valid result structure with NA values
                    return(list(
                        model_summary = list(
                            gene_id = gene,
                            gene_name = gene_annot$gene_name[gene_annot$gencodeID == gene][1],
                            gene_type = gene_annot$gene_type[gene_annot$gencodeID == gene][1],
                            network = network_name,
                            module = module_name,
                            alpha = NA,
                            n_snps_in_window = NA,
                            n_snps_in_model = NA,
                            best_lambda = NA,
                            cv_R2_pc1_avg = NA,
                            cv_R2_pc1_sd = NA,
                            cv_rho_pc1_avg = NA,
                            cv_rho_pc1_se = NA,
                            cv_rho_pc1_avg_squared = NA,
                            cv_zscore_est_pc1 = NA,
                            cv_zscore_pval_pc1 = NA,
                            cv_pval_est_pc1 = NA,
                            cv_R2_gene_avg = NA,
                            cv_R2_gene_sd = NA,
                            cv_rho_gene_avg = NA,
                            cv_rho_gene_se = NA,
                            cv_rho_gene_avg_squared = NA,
                            cv_zscore_est_gene = NA,
                            cv_zscore_pval_gene = NA,
                            cv_pval_est_gene = NA,
                            cor_pc1_all_data_pred = NA,
                            cor_gene_all_data_pred = NA,
                            adj_rsq_gene_all_data = NA,
                            adj_rsq_pc1_all_data = NA,
                            pval_gene_all_data = NA,
                            pval_pc1_all_data = NA,
                            sign_inversion = NA
                        ),
                        weighted_snps_info = NULL
                    ))
                })
            }, mc.cores = config$n_cores, mc.preschedule = FALSE)
            
            names(gene_models) = module$genes
            gene_models        = gene_models[!is.na(gene_models)]
            return(gene_models)
        })
        network_results = network_results[!is.na(network_results)]
       
        
        # Create extra and weights files
        model_summaries = do.call(rbind, lapply(network_results, function(module_results) {
            do.call(rbind, lapply(module_results, function(gene_result) {
                if (!is.null(gene_result) & class(gene_result) == "list") {
                  data.frame(t(gene_result$model_summary), stringsAsFactors = FALSE)
                }
            }))
        }))
        # Convert in numeric columns 
        model_summaries[, 6:(ncol(model_summaries) - 1)] = lapply(model_summaries[, 6:(ncol(model_summaries) - 1)], as.numeric)
        
        ## Compute qvalues (as EpiXcan)
        #model_summaries$qval_gene = qvalue(model_summaries$pval_gene_all_data)$qvalues
        model_summaries$qval_gene  = NA
        
        colnames(model_summaries) = c(
          "gene_id", "gene_name", "gene_type", "network", "module", "alpha",
          "n_snps_in_window", "n_snps_in_model", "best_lambda",
          "cv_R2_pc1_avg", "cv_R2_pc1_sd", "cv_rho_pc1_avg", "cv_rho_pc1_se", "cv_rho_pc1_avg_squared", 
          "cv_zscore_est_pc1", "cv_zscore_pval_pc1", "cv_pval_est_pc1", "cv_R2_gene_avg", "cv_R2_gene_sd",          
          "cv_rho_gene_avg" , "cv_rho_gene_se", "cv_rho_gene_avg_squared", "cv_zscore_est_gene", "cv_zscore_pval_gene",    
          "cv_pval_est_gene", "cor_pc1_all_data_pred", "cor_gene_all_data_pred", "adj_rsq_gene_all_data", "adj_rsq_pc1_all_data",   
          "pval_gene_all_data", "pval_pc1_all_data", "sign_inversion", "qvalue_gene_all_data"                   
        )
        
        # Save model summary
        write.table(
            model_summaries,
            file = file.path(summary_output_dir, paste0(network_name, ".extra.txt")),
            quote = FALSE,
            row.names = FALSE,
            sep = "\t"
          )
        
        # Create weights file
        weights_data = do.call(rbind, lapply(network_results, function(module_results) {
            do.call(rbind, lapply(module_results, function(gene_result) {
                if (!is.null(gene_result) & class(gene_result) == "list") {
                    gene_result$weighted_snps_info
                }
            }))
        }))
        # Save weights
        write.table(
            weights_data,
            file = file.path(weights_output_dir, paste0(network_name, ".weights.txt")),
            quote = FALSE,
            row.names = FALSE,
            sep = "\t"
        )
        
        # Save full results
        #save(network_results,
         #    file = file.path(config$network_dir, 
         #                     paste0(network_name, "_models.RData")))
    }
}

#################################################
## Run Pipeline
#################################################


for(brain_region in CONFIG$brain_regions) {
    print(paste0("Processing brain region: ", brain_region))
    process_brain_region(brain_region, CONFIG)
}

 

