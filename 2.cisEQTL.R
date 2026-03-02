#####################################################################
# Brain Tissue cis-eQTL Analysis Pipeline
#####################################################################
# This script performs cis-eQTL analysis for brain tissue.
# It processes expression and genotype data to:
# 1. Perform cross-validation for eQTL detection
# 2. Rank and prune eQTLs based on significance
# 3. Save results for each chromosome
#####################################################################

################################################
# 1. Setup and Configuration
################################################
# Load required libraries and set up configuration
required_packages = c(
    "data.table", "dplyr", "parallel", "MASS", 
    "sfsmisc", "matrixStats", "purrr"
)

# Install missing packages and load all required ones
for(pkg in required_packages) {
    if(!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
        library(pkg, character.only = TRUE)
    }
}
    
## Set working directory
# Set working directory to project root
# setwd("path/to/project/root")


#################################################
## Helper Functions
#################################################
source("./STEP2.Training_cis/cisEQTL_helper_functions.R")

#####################
# Configuration
#################
base_dir = getwd()

################################################
# 2. Configuration Settings
################################################
CONFIG = list(
    rnaseq              = paste0(base_dir, "/rnaseq/your_rnaseq.RData"),  # Rdata expression, covariates & annotation file
    genotype_by_chr_dir = paste0(base_dir, "/STEP2.Training_CIS/genotype_by_chr/"),  # Directory containing chr*_genotype.rds files
    output_dir          = file.path(base_dir, "/STEP2.Training_CIS/output/"),
    data_dir            = file.path(base_dir, "/data/"),
    n_cores             = 1,  
    cis_window          = 1e6,  # cis-window size in bp (1Mb)
    chromosomes         = 1:22,  # Chromosomes to analyze
    brain_regions       = c("dlpfc"),
    batch_size          = 80  # Number of genes to process in parallel
)

#################################################
## Usage 
#################################################

# Create output directory if it doesn't exist
dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)

# Load expression data
data = get(load(CONFIG$rnaseq))

################################################
# 3. Main Analysis Loop
################################################
# Process each brain region
for(brain_region in CONFIG$brain_regions) {
    message(sprintf("Processing brain region: %s", brain_region))
    
    # Extract and prepare data for current brain region
    expression_data = as.data.frame(data[[brain_region]]$assays$expression)
    covariates      = as.data.frame(data[[brain_region]]$colData)
    gene_annotation = as.data.frame(data[[brain_region]]$rowData)
    
    # Clean chromosome information
    gene_annotation$chr = as.numeric(gsub("chr","",gene_annotation$chr))
    gene_annotation     = gene_annotation[!is.na(gene_annotation$chr),] # Remove chrX, chrY
    
    ## Add TSS
    gene_annotation     = gene_annotation %>% mutate(tss = case_when(strand=="+" ~ start, strand=="-" ~ end))
    
    # Process each chromosome
    for(chr in CONFIG$chromosomes) {
        message(sprintf("Processing chromosome %d", chr))
        
        # Load and validate chromosome data
        chr_file = file.path(CONFIG$genotype_by_chr_dir, brain_region, sprintf("chr%d_genotype.rds", chr))
        if(!file.exists(chr_file)) {
            warning(sprintf("Genotype file not found for chromosome %d, skipping...", chr))
            next
        }
        
        # Get cis-pairs and run analysis
        chr_data = readRDS(chr_file)
        cis_data = filter_cis_pairs(
            gene_annot = gene_annotation,
            snp_annot  = chr_data$snp_info,
            chromosome = chr,
            cis_window = CONFIG$cis_window
        )
        
        if(!is.null(cis_data)) {
            # Process genes in batches
            gene_ids = names(cis_data$pairs)
            n_batches = ceiling(length(gene_ids) / CONFIG$batch_size)
            
            for(batch in 1:n_batches) {
                start_idx   = (batch - 1) * CONFIG$batch_size + 1
                end_idx     = min(batch * CONFIG$batch_size, length(gene_ids))
                batch_genes = gene_ids[start_idx:end_idx]
                
                message(sprintf("Processing batch %d/%d (genes %d-%d)", 
                              batch, n_batches, start_idx, end_idx))
                
                # Process batch in parallel using mclapply
                batch_results = mclapply(batch_genes, function(gene_id) {
                    tryCatch({
                        cis_snps = cis_data$pairs[[gene_id]]
                        if(length(cis_snps) == 0) return(NULL)
                        
                        results = run_rlm_analysis(
                            gene        = gene_id,  
                            expression  = expression_data[,gene_id,drop=FALSE],
                            genotypes   = chr_data$genotype[,cis_snps],
                            snp_annot   = chr_data$snp_info,
                            covariates  = covariates,
                            n_cores     = 1,  # Use 1 core per gene since we're parallelizing at batch level
                            output_dir  = CONFIG$output_dir,
                            data_dir    = CONFIG$data_dir,
                            cv_fold_id  = "your_fold_id",
                            brain.region = brain_region
                        )
                        
                        if(!is.null(results)) {
                            results$gene       = gene_id
                            results$snp        = rownames(results)
                            results$chromosome = chr
                        }
                        return(results)
                    }, error = function(e) {
                        message(sprintf("Error processing gene %s: %s", gene_id, e$message))
                        return(NULL)
                    })
                }, mc.cores = CONFIG$n_cores, mc.preschedule = TRUE)
                
                # Combine batch results
                batch_results = do.call(rbind, batch_results[!sapply(batch_results, is.null)])
                
                # Save intermediate results
                if(!is.null(batch_results) && nrow(batch_results) > 0) {
                    output_file = file.path(CONFIG$output_dir, brain_region,
                                          sprintf("chr%d_batch%d_significant_ciseQTL.Rds", chr, batch))
                    saveRDS(batch_results, output_file)
                }
                
                # Clear memory
                rm(batch_results)
                gc()
            }
            
            # Combine all batch results
            all_results = list.files(
                path = file.path(CONFIG$output_dir, brain_region),
                pattern = sprintf("chr%d_batch.*_significant_ciseQTL.Rds", chr),
                full.names = TRUE
            )
            
            if(length(all_results) > 0) {
                chr_results = do.call(rbind, lapply(all_results, readRDS))
                
                # Save final results
                output_file = file.path(CONFIG$output_dir, brain_region,
                                      sprintf("chr%d_significant_ciseQTL.Rds", chr))
                saveRDS(chr_results, output_file)
                
                # Clean up intermediate files
                unlink(all_results)
            }
        }
        
        rm(chr_data, cis_data)
        gc()
    }
   
    message(sprintf("Region %s completed", brain_region))
    gc()
}
