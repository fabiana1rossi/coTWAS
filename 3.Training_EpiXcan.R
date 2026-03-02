#####################################################################
# Brain Tissue Gene Expression Prediction Training Pipeline
#####################################################################

# Required packages
required_packages = c(
  "data.table", "dplyr", "glmnet", "parallel",
  "limma"
)

# Install and load required packages
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


## Set working directory
# Set working directory to project root
# setwd("")


###################
## HELPER FUNCTIONS
###############
source("./STEP2.Training_cis/training_helper_functions_EpiXcan.R")


# Configuration settings
base_dir = getwd()

CONFIG   = list(
  # Input/Output directories
  rnaseq              = paste0(base_dir, "/rnaseq/your_rnaseq.RData"),  # Rdata expression, covariates & annotation file
  genotype_by_chr_dir = paste0(base_dir, "/STEP2.Training_cis/genotype_by_chr_LIBD/"),  # Directory containing chr*_genotype.rds files
  output_dir          = file.path(base_dir, "/STEP2.Training_cis/output_EpiXcan/"),
  eqtls_dir           = file.path(base_dir, "/STEP2.Training_cis/output/"),
  prior_dir           =  "../../qtlBHM/results_LIBD_annotations_qtlBHM/", ## your qtlBHM results folder
  annot_dir           =  "../../qtlBHM/annot/", ## your qtlBHM annotation folder
  data_dir            = file.path(base_dir, "/data/"),

  # Data parameters
  brain_regions = c("dlpfc"),
  chromosomes   = 1:22,

  # Model parameters
  n_folds    = 4,  # Number of CV folds
  cis_window = 1e6,  # Size of cis window in base pairs
  alpha      = 0.5,  # Elastic net mixing parameter

  # Computation parameters
  n_cores = 20,
  seed    = 070892,

  # Covariate parameters
  covariate_selection = list(
    demographic = c("Age", "Sex", "Dx"),
    technical   = c("RIN", "mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate"),
    batch       = c("Protocol", "Dataset"),
    pcs         = paste0("PC", 1:5)
  ),

  # Output parameters - Updated to match actual output from helper functions
  model_summary_cols = c(
    'gene_id', 'gene_name', 'gene_type', 'chromosome',
    'alpha', 'n_snps_in_window', 'n_snps_in_model',
    'best_lambda', 'inner_R2_avg', 'inner_R2_sd', 'inner_pval_est',
    'inner_rho_avg', 'inner_rho_se', 'inner_rho_zscore', 'inner_rho_avg_squared',
    'inner_zscore_pval', 'cv_R2_gene_avg', 'cv_R2_gene_sd',
    'cv_rho_gene_avg', 'training_R2', 'cv_rho_gene_se', 'cv_rho_gene_avg_squared',
    'cv_zscore_est_gene', 'cv_zscore_pval_gene', 'cv_pval_est_gene',
    'cor_gene_all_data_pred', 'adj_rsq_gene', 'rmse_avg',
    'pval_gene', 'cv_adj_r2_gene_avg', 'cv_adj_r2_gene_sd', 'cv_pval_lm_gene_avg'
  ),

  weights_cols = c(
    'gene', 'gene_name', 'gene_start', 'gene_end', 'strand', 'chr', 'pos',
    'distTSS', 'rsid', 'Annotation', 'Prior', 'penalty_factor',
    'refAllele', 'effectAllele', 'weight', 'alpha'
  ),

  # File patterns
  file_patterns = list(
    gene_annot    = "%s.gene_annot.txt",
    expression    = "%s.ranknorm_expression.txt",
    covariates    = "%s_combined.covariates.txt",
    geno          = "chr%d_genotype.Rds",
    eqtl          = "chr%d_signifcant_ciseQTL.Rds",
    prior         = "results_%s_variant_posterior.txt",
    weights       = "%s_chr%d_weights.txt",
    summary       = "%s_chr%d_model_summaries.txt",
    model_summary = "%s_chr%d_model_summary.txt",
    cv_folds      = "CIS_cv_4fold_%s_ids.RData"
  )
)


#' Main pipeline function
run_training_pipeline = function(CONFIG) {
  # Load data
  data = get(load((CONFIG$rnaseq)))
  
  for(brain_region in CONFIG$brain_regions) {
    message(sprintf("Processing brain region: %s", brain_region))
    
    # Create output directories
    if(!file.exists(file.path(CONFIG$output_dir, "summary", brain_region))){
      dir.create(file.path(CONFIG$output_dir, "summary", brain_region), recursive=TRUE)}
    if(!file.exists(file.path(CONFIG$output_dir, "weights", brain_region))){
      dir.create(file.path(CONFIG$output_dir, "weights", brain_region), recursive=TRUE)}
    
    ## Get brain region data
    data_region = data[[brain_region]]
    
    ## Open qtlBHM annotation results
    results_annotations = read.delim(file.path(CONFIG$prior_dir, sprintf("results_%s_annotations.txt", toupper(brain_region))))
    
    ## Open annotation BED
    annotation_bed = read.delim(file.path(CONFIG$annot_dir, sprintf("%s_25_imputed12marks_mnemonics.bed", brain_region)), header=FALSE)
    # Set column names for annotation file (BED format)
    colnames(annotation_bed) = c("Chr", "Start", "End", "Annotation")
    
    # Clean chromosome names to match snp_map format
    annotation_bed$Chr = gsub("chr", "", annotation_bed$Chr)
    
    # Create genomic position key for annotation
    annotation_bed$pos_key = paste0(annotation_bed$Chr, ":", annotation_bed$Start)
    
    # Process each chromosome
    for(chr in CONFIG$chromosomes) {
      process_chromosome(brain_region, chr, data_region, CONFIG, annot_results = results_annotations, annot_bed = annotation_bed)
    }
    
    # Combine results across chromosomes
    message("Combining results across chromosomes...")
    
    # Combine summary files
    all_summaries = do.call(rbind, lapply(CONFIG$chromosomes, function(chr) {
      summary_file = get_file_path("summary", brain_region, chr, CONFIG)
      if(file.exists(summary_file) && file.info(summary_file)$size > 2) { ## If file exists and it is not empty
        read.delim(summary_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
      }
    }))
    
    # Save combined summary
    write.table(
      all_summaries,
      file = file.path(CONFIG$output_dir, "summary",brain_region,
                       paste0(brain_region, "_all_chromosomes_summary.txt")),
      quote = FALSE,
      row.names = FALSE,
      sep = "\t"
    )
    
    # Combine weights files
    all_weights = do.call(rbind, lapply(CONFIG$chromosomes, function(chr) {
      weights_file = get_file_path("weights", brain_region, chr, CONFIG)
      if(file.exists(weights_file) && file.info(weights_file)$size > 2) { ## If file exists and it is not empty
        read.delim(weights_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
      }
    }))
    
    # Save combined weights
    write.table(
      all_weights,
      file = file.path(CONFIG$output_dir, "weights", brain_region,
                       paste0(brain_region, "_all_chromosomes_weights.txt")),
      quote = FALSE,
      row.names = FALSE,
      sep = "\t"
    )
    
    # Optionally, remove individual chromosome files
    if(TRUE) { # Set to TRUE if you want to clean up individual files
      for(chr in CONFIG$chromosomes) {
        unlink(get_file_path("summary", brain_region, chr, CONFIG))
        unlink(get_file_path("weights", brain_region, chr, CONFIG))
      }
    }
  }
}

################
# Run pipeline 
#################
run_training_pipeline(CONFIG)
