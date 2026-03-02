# This script identifies SNPs that are present in the EpiXcan database but missing from the target genotype data
# It is typically used as a preprocessing step for genetic analysis to identify variants that need proxy SNPs

# Load required packages
if (!require("data.table")) install.packages("data.table")
if (!require("RSQLite")) install.packages("RSQLite") 
if (!require("dplyr")) install.packages("dplyr")
if (!require("limma")) {
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("limma")
}

library(data.table)
library(RSQLite)
library(dplyr) 
library(limma)

# Set working directory to project root
# setwd("path/to/project/root")  ## your path here

## Plink path 
plink = "../../tools/plink"


#####

# Configuration
base_dir = getwd()
CONFIG = list(
  out_dir   = file.path(base_dir, "STEP4.LDproxy_Testing_Geno/LD_output"),
  db_dir    = file.path(base_dir, "STEP3.Training_MODULE/output/database/"),
  geno_dir  = file.path(base_dir, "genotype"),
  geno_name = "your_geno_prefix"  
) 


# Load target SNPs from the PLINK binary format (.bim file)
# The .bim file contains the genetic variants present in your dataset
load_target_snps = function(bim_file) {
  fread(bim_file)
}

# Identify SNPs that are in the reference (EpiXcan) but not in the target dataset
# These SNPs will need proxy variants for downstream analysis
find_mismatches = function(ref_snps, target_snps) {
  setdiff(ref_snps$ID, target_snps$V2)
}

# Main execution function
main = function() {
  start_time = Sys.time()
  
  # Step 1: Load SNPs from  database (reference)
  ref = get(load(file.path(CONFIG$db_dir, "all_weights_LD.RData")))
  
  # Step 2: Load SNPs from target genotype data
  target = load_target_snps(file.path(CONFIG$geno_dir, paste0(CONFIG$geno_name, ".bim")))
  
  # Step 3: Identify SNPs present in EpiXcan but missing in target data
  mismatch = find_mismatches(ref, target)
  
  # Step 4: Save mismatched SNPs to a file for further processing
  # These SNPs will need proxy variants in subsequent analysis steps
  write.table(mismatch, 
              file = {
                if (!dir.exists(CONFIG$out_dir)) dir.create(CONFIG$out_dir)
                file.path(CONFIG$out_dir, paste0(CONFIG$geno_name, "_mismatch.txt"))
              },
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  print(Sys.time() - start_time)
}

main()


