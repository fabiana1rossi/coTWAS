
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
  db_dir1    = file.path(base_dir, "STEP2.Training_CIS/output/database/"),
  db_dir2    = file.path(base_dir, "STEP2.Training_CIS/output_EpiXcan/database/"),
  geno_dir  = file.path(base_dir, "genotype"),
  geno_name = "your_genotype_name" # GTeX.v9.phased.geno.maf.bfile.updated.EUR   CMC_EUR.IBD.geno0.02.hwe10e6.mind0.02.maf0.01.hg20  LIBD_TopMed.geno.maf.EUR.chr1_22.noindels.overlap.GTEX
) 


# Load target SNPs from the PLINK binary format (.bim file)
# The .bim file contains the genetic variants present in your dataset
load_target_snps = function(bim_file) {
  fread(bim_file)
}

# Identify SNPs that are in the reference (EpiXcan) but not in the target dataset
# These SNPs will need proxy variants for downstream analysis
find_mismatches = function(ref_snps, target_snps) {
  setdiff(ref_snps$ref, target_snps$V2)
}

# Main execution function
main = function() {
  start_time = Sys.time()
  
  # Step 1: Load SNPs from  database (reference)
  ref = readRDS(file.path(CONFIG$db_dir1, "weights_LD.rds"))   ## dataframe containing gene-snps-ref_allele-alt_allele-POS-CHR  (from SQLite DBs)
  ref = distinct(ref)
  
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


