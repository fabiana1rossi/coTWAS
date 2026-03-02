#' CIS Genotype Matching Script --> Output will be used in INGENE training --> imputed coexpression partners values
#' 
#' This script performs genotype matching for CIS analysis by:
#' 1. Loading and processing genotype data from bim files
#' 2. Matching SNPs with CIS reference database
#' 3. Handling allele flipping and complementing
#' 4. Converting matched genotypes to VCF format
#' 
#' Required libraries: data.table, magrittr, purrr, RSQLite
#' Required external tools: plink2

# Load required libraries
require(data.table)
require(magrittr)
require(purrr)
require(RSQLite)


## set working directory
# Set working directory to project root
# setwd("path/to/project/root")

# Get base directory for relative paths
base_dir = getwd()

# Configuration settings
CONFIG    = list(
  inputs = list(
    file       = "training_genotype", 
    directory  = "./genotype/",
    target     = "training_genotype",
    regions    = c("dlpfc")
  ),
  output  = list(
    directory  = "./genotype/CIS.matched.geno/"
  ),
  tools   = list(
    plink2_path = "../../tools/plink2"
  ),
  database  = list(
    base_path = "./STEP2.Training_CIS/output/database/",
    suffix    = "CIS_%s_testing_genotype_name"
  ),
  temp_files = list(
    directory        = "./temp/",
    to_extract       = "CIS.to.extract.txt",
    to_update_allele = "CIS.to.update.allele.txt",
    to_swap_allele   = "CIS.to.swap.allele.txt"
  )
)

#' Nucleotide complement function
#' Returns the complementary nucleotide for DNA base pairs
complement = function(x){
  switch (x,
          "A" = "T",
          "C" = "G",
          "T" = "A",
          "G" = "C",
          return(NA)
  )
} 

#' Main processing function for CIS matching
#' @param input Base name of input file
#' @param input_dir Directory containing input files
#' @param target Target dataset identifier
#' @param output_dir Directory for output files
#' @param config Configuration list
process_cis_match = function(brain_region, input, input_dir, target, output_dir,  config) {
  
  print(sprintf("Processing input: %s in %s ", input, brain_region))
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    print(paste("Created output directory:", output_dir))
  }
  
  # STEP 1: Load and filter bim file
  bim = fread(paste0(input_dir, input, '.bim')) %>% 
    setnames(.,colnames(.)[c(1,2,4,5,6)],c("CHR","ID","POS","ALT",'REF'))
  
  # Remove ambiguous SNPs
  condition = (bim$REF == "A" & bim$ALT == "T") | 
             (bim$REF == "T" & bim$ALT == "A") |
             (bim$REF == "G" & bim$ALT == "C") | 
             (bim$REF == "C" & bim$ALT == "G")
  bim = bim[!condition,]
  bim$CHR = as.numeric(bim$CHR)
  
  # STEP 2: Load CIS reference database
  print("matching with CIS ref")
  sqlite.driver = dbDriver("SQLite")
  db_path = file.path(config$database$base_path,brain_region, paste0(sprintf(config$database$suffix,brain_region), ".db"))
  db      = dbConnect(sqlite.driver, dbname = db_path)
  
  df     = dbReadTable(db,"weights")
  df$ID  = df$rsid
  df$CHR = df$chr
  df$REF = df$ref_allele
  df$EFF = df$eff_allele
  df     = df[,c("ID","CHR","REF","EFF")]
  dbDisconnect(db)
  
  # STEP 3: Match SNPs
  info = merge(bim, df, by=c("ID")) %>% .[!duplicated(ID)]
  
  # Identify different types of SNP matches
  info.match  = info[REF.x == REF.y & ALT == EFF, ID]  # Direct matches
  com.snps    = info[sapply(REF.x, complement) == REF.y &
                   sapply(ALT, complement) == EFF, ID]  # Complement matches
  recode.snps = info[REF.x==EFF & ALT==REF.y, ID]     # Needs recoding
  com.recode  = info[sapply(REF.x, complement) == EFF &
                    sapply(ALT, complement) == REF.y, ID] # Needs complement & recoding
  
  # Handle mismatching SNPs
  mismatch = info$ID[!(info$ID %in% info.match |
                        info$ID %in% com.snps | 
                        info$ID %in% recode.snps |
                        info$ID %in% com.recode)]
  print(paste0(length(mismatch)," SNPs unresolved"))
  
  # Remove mismatching SNPs
  info = info[!info$ID %in% mismatch,c("ID","REF.x","ALT","REF.y","EFF")]
  
  # Prepare files for PLINK
  to.flip          = info[info$ID %in% com.snps,]
  to.flip.recode   = info[info$ID %in% com.recode,c("ID","REF.x","ALT","EFF","REF.y")]
  colnames(to.flip.recode) = colnames(to.flip)
  to.update.allele = rbind(to.flip,to.flip.recode)
  to.swap          = info[info$ID %in% c(recode.snps,com.recode),]
  to.extract       = info$ID
  
  # Create temp directory if needed
  if (!dir.exists(config$temp_files$directory)) {
    dir.create(config$temp_files$directory, recursive = TRUE)
  }
  
  # Set up file paths
  extract_path = file.path(config$temp_files$directory, config$temp_files$to_extract)
  update_path  = file.path(config$temp_files$directory, config$temp_files$to_update_allele)
  swap_path    = file.path(config$temp_files$directory, config$temp_files$to_swap_allele)
  
  # Write extract file
  cat(to.extract, sep = "\n", file = extract_path)
  
  # STEP 4: Execute PLINK commands
  if(!length(Reduce(union,list(com.snps,recode.snps,com.recode)))){
    cmd.1 = paste0(config$tools$plink2_path, 
                  " --bfile ", input_dir, input,
                  " --extract ", extract_path,
                  " --export vcf --out ", output_dir, input, "_",brain_region,"_CIS.matched")
    system(cmd.1)
    print("extraction and vcf conversion completed")
  }
  
  if(length(union(com.snps,com.recode)) & !length(union(recode.snps,com.recode))){
    write.table(to.update.allele, sep = "\t", row.names = F, col.names = F, quote = F,
                file = update_path)
    cmd.1 = paste0(config$tools$plink2_path, " --bfile ", input_dir, input,
                  " --extract ", extract_path, " --make-bed --out ", output_dir, input, ".extract")
    cmd.2 = paste0(config$tools$plink2_path, " --bfile ", output_dir, input, ".extract",
                  " --update-alleles ", update_path, " --recode vcf --out ", output_dir, input, "_",brain_region,"_CIS.matched")
    system(cmd.1)
    print("extraction completed")
    system(cmd.2)
    print("flipping and vcf conversion completed")
  }
  
  if(!length(union(com.snps,com.recode)) & length(union(recode.snps,com.recode))){
    write.table(to.swap, sep = "\t", row.names = F, col.names = F, quote = F,
                file = swap_path)
    cmd.1 = paste0(config$tools$plink2_path, " --bfile ", input_dir, input,
                  " --extract ", extract_path, " --make-bed --out ", output_dir, input, ".extract")
    cmd.2 = paste0(config$tools$plink2_path, " --bfile ", output_dir, input, ".extract",
                  " --alt1-allele ", swap_path, " 5 1 --recode vcf --out ", output_dir, input, "_",brain_region,"_CIS.matched")
    system(cmd.1)
    print("extraction completed")
    system(cmd.2)
    print("swapping and vcf conversion completed")
  }
  
  if(length(union(com.snps,com.recode)) & length(union(recode.snps,com.recode))){
    write.table(to.update.allele, sep = "\t", row.names = F, col.names = F, quote = F,
                file = update_path)
    write.table(to.swap, sep = "\t", row.names = F, col.names = F, quote = F,
                file = swap_path)
    cmd.1 = paste0(config$tools$plink2_path, " --bfile ", input_dir, input,
                  " --extract ", extract_path, " --make-bed --out ", output_dir, input, ".extract")
    cmd.2 = paste0(config$tools$plink2_path, " --bfile ", output_dir, input, ".extract",
                  " --update-alleles ", update_path, " --make-bed --out ", output_dir, input, ".update.allele")
    cmd.3 = paste0(config$tools$plink2_path, " --bfile ", output_dir, input, ".update.allele",
                  " --alt1-allele ", swap_path, " 5 1 --recode vcf --out ", output_dir, input, "_",brain_region,"_CIS.matched")
    system(cmd.1)
    system(cmd.2)
    print("flipping completed")
    system(cmd.3)
    print("swapping and vcf conversion completed")
  }
  

}

#' Create required directories
create_directories = function(config) {
  dirs_to_create = c(
    config$inputs$directory,
    config$output$directory,
    config$temp_files$directory,
    dirname(config$database$base_path)
  )
  
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      print(paste("Created directory:", dir))
    }
  }
}

# Initialize directories
create_directories(CONFIG)

# Execute main processing
for(brain_region in CONFIG$inputs$regions){
  process_cis_match(
    brain_region = brain_region,
    input        = CONFIG$inputs$file,
    input_dir    = CONFIG$inputs$directory,
    target       = CONFIG$inputs$target,
    output_dir   = CONFIG$output$directory,
    config       = CONFIG
  )
  
  
  
  # Clean up temporary files
  if (dir.exists(CONFIG$temp_files$directory)) {
    unlink(list.files(CONFIG$temp_files$directory, full.names = TRUE))
  }
  # Clean up temporary files
  if (dir.exists(CONFIG$output$directory)) {
    unlink(list.files(CONFIG$output$directory, full.names = TRUE)[!grepl("vcf",list.files(CONFIG$output$directory, full.names = TRUE))])
  }
  
  
  
}
