#####################################################################
# Create SNP Annotation Lookup Table 
# This script processes all SNPs once to create a lookup table
# that maps SNPs to their annotations and prior weights
# OPTIMIZATION: Processes SNPs in batches to avoid memory issues
#####################################################################

# Required packages
required_packages = c(
  "data.table", "dplyr", "parallel"
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
# setwd("path/to/project/root")

# Configuration settings
base_dir = getwd()

CONFIG = list(
  # Input directories
  genotype_by_chr_dir = paste0(base_dir, "/STEP2.Training_CIS/genotype_by_chr_GTEx/"),
  prior_dir           = "../../qtlBHM/results_annotations_qtlBHM/", ## your qtlBHM results folder
  annot_dir           = "../../qtlBHM/annot/", ## folder containing RoadMap Annotations
  output_dir          = paste0(base_dir, "/data/"),
  
  # Data parameters
  brain_regions = c("dlpfc"),
  chromosomes   = 1:22,
  
  # Memory optimization parameters
  batch_size = 1000,  # Process SNPs in batches of 1000
  n_cores = 10,       # Reduced number of cores to avoid memory overload
  
  # Computation parameters
  gc_frequency = 5    # Run garbage collection every 5 batches
)

#' Helper function to get file paths
get_file_path = function(file_type, brain_region=NULL, chr=NULL, CONFIG) {
  pattern = CONFIG$file_patterns[[file_type]]
  
  if(file_type %in% c("geno")) {
    return(sprintf(file.path(CONFIG$genotype_by_chr_dir, brain_region,pattern), chr))
  }
}

# Add file patterns to CONFIG
CONFIG$file_patterns = list(
  geno = "chr%d_genotype.Rds"
)

#' Process a batch of SNPs to create annotation lookup
#' @param snp_batch Batch of SNP indices to process
#' @param snp_annot SNP annotation data
#' @param annot_bed Annotation BED data
#' @param annot_results Annotation results data
#' @return Data frame with SNP annotations for the batch
process_snp_batch = function(snp_batch, snp_annot, annot_bed, annot_results) {
  batch_results = list()
  
  for(i in seq_along(snp_batch)) {
    idx = snp_batch[i]
    snp_chr = snp_annot$Chromosome[idx]
    snp_pos = snp_annot$Position[idx]
    
    # Find matching annotation regions
    matching_annot = annot_bed[as.numeric(annot_bed$Chr) == snp_chr & 
                                 as.numeric(annot_bed$Start) <= snp_pos & 
                                 as.numeric(annot_bed$End) >= snp_pos, ]
    matching_annot = na.omit(matching_annot)
    
    if(nrow(matching_annot) > 0) {
      # If multiple annotations, take the first one
      annot_value = matching_annot$Annotation[1]
    } else {
      annot_value = NA
    }
    
    # Find prior from results_annotations
    # If there is annotation, assign prior=weight; otherwise, assign prior=0
    snp_id = snp_annot$SNP[idx]
    
    if(!is.na(annot_value)) {
      prior_row = annot_results[annot_results$Annotation == annot_value, ]
      prior_row = na.omit(prior_row)
      
      if(nrow(prior_row) > 0) {
        if(nrow(prior_row) > 1) {
          ## Select the row with the highest weight
          prior_row = prior_row[which.max(prior_row$Weight), ]
        }
        prior_value = prior_row$Weight
      } else {
        prior_value = 0  # No matching annotation found in results
      }
    } else {
      prior_value = 0  # No annotation found for this SNP
    }
    
    # Create row for output
    output_row = data.frame(
      Chr = snp_chr,
      Pos = snp_pos,
      Prior = prior_value,
      RSID_dbSNP137 = snp_id,
      Annotation = annot_value,
      stringsAsFactors = FALSE
    )
    
    batch_results[[i]] = output_row
  }
  
  # Combine batch results
  if(length(batch_results) > 0) {
    return(do.call(rbind, batch_results))
  } else {
    return(data.frame())
  }
}


#' Process a single chromosome to create SNP annotation lookup
#' @param chr Chromosome number
#' @param brain_region Brain region name
#' @param CONFIG Configuration parameters
process_chromosome_annotations = function(chr, brain_region, CONFIG) {
  message(sprintf("Processing chromosome %d for %s", chr, brain_region))
  
  # Output file path
  output_file = file.path(CONFIG$prior_dir, sprintf("%s_chr%s_annotations.Rds", brain_region, chr))
  
  # Skip if results exist
  if(file.exists(output_file)) {
    message("Results already exist, skipping...")
    return(NULL)
  }
  
  # Load chromosome-specific genotype data
  geno_file = get_file_path("geno", chr=chr, brain_region = brain_region, CONFIG = CONFIG)
  if(!file.exists(geno_file)) {
    warning(sprintf("Genotype file not found for chromosome %d: %s", chr, geno_file))
    return(NULL)
  }
  
  message("Loading genotype data...")
  genotypes_snp_annot = readRDS(geno_file)
  if(is.null(genotypes_snp_annot) || !all(c("genotype", "snp_info") %in% names(genotypes_snp_annot))) {
    stop(sprintf("Invalid genotype data format for chromosome %d", chr))
  }
  
  snp_annot = genotypes_snp_annot$snp_info
  
  # Clear genotype data from memory to save space
  rm(genotypes_snp_annot)
  gc()
  
  # Load annotation results
  message("Loading annotation results...")
  ## dACC and sACC have the same annotation
  results_file = ifelse(brain_region %in% c("sACC","dACC"), file.path(CONFIG$prior_dir, sprintf("results_%s_annotations.txt", brain_region)),
                                                                      file.path(CONFIG$prior_dir, sprintf("results_%s_annotations.txt", toupper(brain_region))))
  if(!file.exists(results_file)) {
    stop(sprintf("Annotation results file not found: %s", results_file))
  }
  annot_results = read.delim(results_file)
  
  # Load annotation BED file
  message("Loading annotation BED file...")
  bed_file = file.path(CONFIG$annot_dir, sprintf("%s_25_imputed12marks_mnemonics.bed", brain_region))
  if(!file.exists(bed_file)) {
    stop(sprintf("Annotation BED file not found: %s", bed_file))
  }
  annot_bed = read.delim(bed_file, header=FALSE)
  colnames(annot_bed) = c("Chr", "Start", "End", "Annotation")
  
  # Clean chromosome names
  annot_bed$Chr = gsub("chr", "", annot_bed$Chr)
  
  message(sprintf("Processing %d SNPs for chromosome %d in batches of %d", 
                 nrow(snp_annot), chr, CONFIG$batch_size))
  
  # Process SNPs in batches
  total_snps = nrow(snp_annot)
  n_batches = ceiling(total_snps / CONFIG$batch_size)
  all_results = list()
  
  for(batch in 1:n_batches) {
    start_idx = (batch - 1) * CONFIG$batch_size + 1
    end_idx = min(batch * CONFIG$batch_size, total_snps)
    
    message(sprintf("Processing batch %d/%d (SNPs %d-%d)", 
                   batch, n_batches, start_idx, end_idx))
    
    # Create batch indices
    batch_indices = start_idx:end_idx
    
    # Process batch in parallel
    batch_size_small = ceiling(length(batch_indices) / CONFIG$n_cores)
    batch_chunks = split(batch_indices, ceiling(seq_along(batch_indices) / batch_size_small))
    
    batch_results = mclapply(batch_chunks, function(chunk) {
      process_snp_batch(snp_batch = chunk, snp_annot, annot_bed, annot_results)
    }, mc.cores = CONFIG$n_cores, mc.preschedule = FALSE)
    
    # Combine batch results
    batch_combined = do.call(rbind, batch_results)
    all_results[[batch]] = batch_combined
    
    # Clear batch results from memory
    rm(batch_results, batch_combined)
    
    # Run garbage collection periodically
    if(batch %% CONFIG$gc_frequency == 0) {
      gc()
      message(sprintf("Garbage collection completed after batch %d", batch))
    }
  }
  
  # Combine all results
  message("Combining all batch results...")
  snp_prior = do.call(rbind, all_results)
  
  # Clear intermediate results from memory
  rm(all_results)
  gc()
  
  # Save results
  saveRDS(snp_prior, file = output_file)
  message(sprintf("Saved %d SNP annotations for chromosome %d", nrow(snp_prior), chr))
  
  return(snp_prior)
}

#' Main function to create SNP annotation lookup for all chromosomes
create_snp_annotation_lookup = function(CONFIG) {
  for(brain_region in CONFIG$brain_regions) {
    message(sprintf("Processing brain region: %s", brain_region))
    
    # Process each chromosome
    for(chr in CONFIG$chromosomes) {
      process_chromosome_annotations(chr, brain_region, CONFIG)
      
      # Run garbage collection between chromosomes
      gc()
      message(sprintf("Memory cleanup completed for chromosome %d", chr))
    }
    
    message(sprintf("Completed processing for brain region: %s", brain_region))
  }
}

################
# Run the script
#################
create_snp_annotation_lookup(CONFIG) 
