
#################################################
## Script to prepare genotype data and map SNPs to gene coexpression modules
##
## This script performs:
## 1. Maps SNPs to gene co-expression modules:
##    - Loads gene coexpression module from file and SNP annotation data
##    - For each brain region:
##      * Gets gene sets from expression data
##      * Identifies SNPs within 100kb of genes in modules 
##      * Filters SNPs present in genotype data
##      * Creates region-specific SNP lists for each module in coexpression network
##
## Input files required:
##  - PLINK format genotype files (.bed/.bim/.fam)
##  - Gene coexpression module file (RDS file)
##  - SNP annotation data from 1000G (100kb window, RData file)
##  - RNA-seq expression data (RData file)
##
## Output:
##  - RData files containing SNP lists mapped to modules
##  - Returns a list of SNP mappings for each brain region
#################################################

## Required packages
if (!require("BiocManager")) install.packages("BiocManager")

# CRAN packages
if (!require("dplyr")) install.packages("dplyr")
if (!require("purrr")) install.packages("purrr") 
if (!require("parallel")) install.packages("parallel")
if (!require("limma")) install.packages("limma")

# Bioconductor packages
if (!require("SummarizedExperiment")) BiocManager::install("SummarizedExperiment")
if (!require("SNPRelate")) BiocManager::install("SNPRelate")

library(dplyr)
library(purrr)
library(parallel)
library(limma)
library(SummarizedExperiment)
library(SNPRelate)

## set working directory
# setwd("")


#################################################
## Helper Functions  ##
#################################################
source("./STEP3.Training_MODULE/prepare_mapping_helper_functions.R")
"%&%" = function(a,b) paste(a,b, sep='')



###################
## Configuration ##
###################

# Set paths
base_dir      = getwd()
#plink_version = paste0(base_dir,"/tools/plink2")
geno_dir      = paste0(base_dir, "/genotype/")
rnaseq_dir    = paste0(base_dir, "/rnaseq/")
data_dir      = paste0(base_dir, "/data/")

# Input files
geno_file                  = paste0(geno_dir, "your_genotype_prefix")
geno_name                  = "your_genotype_prefix"
networks_file              = paste0(data_dir, "coexpression_network.rds") ## your coexpression networks
snplist_file               = paste0(data_dir, "100kbp_SNPList.1000G.EUR.RData")
expr_file                  = paste0(rnaseq_dir, "your_rnaseq.RData")

# Output directory
output_dir = paste0(base_dir, "/STEP3.Training_MODULE/output/")

# Parameters
brain_regions = c("dlpfc") ## Specify all available regions 
n_cores       = 10          ## Specify more cores to run code in parallel on ubuntu terminal


#####################
## Main Execution ##
#####################

main = function() {
  
  # Process SNPs and create mapping
  snpList_modules = process_snps(networks_file, snplist_file, expr_file, geno_dir,geno_name, output_dir, brain_regions, n_cores)
  
  # Return results
  return(snpList_modules)
}

# Run pipeline
results = main() 
