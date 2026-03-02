#####################################################################
# Genotype Data Processing Pipeline - Optimized for cis-eQTL analysis
# Brain Region Specific Version
#####################################################################
# This script loads and processes genotype data from PLINK format,
# converts it to GDS, and creates chromosome-specific RDS files for each brain region
# Only cis SNPs are considered for each brain region
#####################################################################

################################################
## Setup and Configuration
#################################################

required_packages = c("SNPRelate", "gdsfmt", "data.table", "dplyr", "bigmemory", "SummarizedExperiment")

for(pkg in required_packages) {
    if(!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
        library(pkg, character.only = TRUE)
    }
}

## Set working directory
# Set working directory to project root
# setwd("")

# Define brain regions
BRAIN_REGIONS = c("dlpfc")

CONFIG = list(
    rnaseq        = file.path(getwd(), "/rnaseq/your_rnaseq.RData"),
    genotype_path = paste0(getwd(), "/genotype/your_genotype_prefix"),
    output_dir    = file.path(getwd(), "STEP2.Training_CIS/genotype_by_chr/"),
    cis_window    = 1e6,  # 1Mb cis window
    brain_regions = BRAIN_REGIONS
)

#################################################
## Helper Functions
#################################################

#' Get genes and their positions from expression data for a specific brain region
#' @param expression_file Path to training expression data
#' @param region Brain region to process
#' @return Data frame with gene information
get_gene_positions = function(expression_file, region) {
    # Load expression data
    expr_data = get(load(expression_file))
    
    # Get region-specific data
    if (!region %in% names(expr_data)) {
        stop(sprintf("Region %s not found in expression data", region))
    }
    
    # Extract gene annotation for the specific region
    gene_info = as.data.frame(expr_data[[region]]$rowData)
    
    # Create data frame with required information
    gene_positions = data.frame(
        gene_id = rownames(gene_info),
        chr     = as.numeric(gsub("chr", "", gene_info$chr)),
        start   = gene_info$start,
        end     = gene_info$end,
        strand  = gene_info$strand,
        stringsAsFactors = FALSE
    )
    
    gene_positions = gene_positions[!is.na(gene_positions$chr),]
    
    ## Add TSS column
    gene_positions = gene_positions %>% mutate(tss = case_when(strand=="+" ~ start, strand=="-" ~ end))
    
    message(sprintf("Processing %d genes for %s", nrow(gene_positions), region))
    return(gene_positions)
}

#' Find SNPs within cis windows of genes
#' @param snp_info SNP information
#' @param gene_positions Gene position information
#' @param window_size Size of cis window in bp
#' @return Vector of SNP IDs within cis windows
find_cis_snps = function(snp_info, gene_positions, window_size) {
    # Initialize empty list to store cis SNPs
    cis_snps_list = list()
    
    # Process each chromosome
    for (chr in unique(gene_positions$chr)) {
        # Get genes and SNPs for this chromosome
        gene_pos_chr = gene_positions[gene_positions$chr == chr,]
        snp_pos_chr  = snp_info[snp_info$chromosome == chr,]
        
        if (nrow(gene_pos_chr) == 0 || nrow(snp_pos_chr) == 0) next
        
        # Create vectors of tss for all genes on chromosome
        gene_tss = gene_pos_chr$tss
        
        # For each SNP, check if it falls within any gene's window
        snp_positions = snp_pos_chr$position
        
        # Find SNPs in cis windows
        in_cis = sapply(snp_positions, function(pos) {
            any((pos >= (gene_tss - 1e6)) & (pos <= (gene_tss + 1e6)))
        })
        
        # Store cis SNPs for this chromosome
        if (any(in_cis)) {
            cis_snps_list[[as.character(chr)]] = snp_pos_chr$snp.id[in_cis]
        }
    }
    
    # Combine all cis SNPs and remove duplicates
    cis_snps = unique(unlist(cis_snps_list))
    
    message(sprintf("Found %d unique cis SNPs across all chromosomes", length(cis_snps)))
    return(cis_snps)
}

#' Process genotype data
#' @param genofile GDS file handle
#' @param snp_ids SNP IDs to process
#' @param sample_ids Sample IDs
#' @return Genotype matrix
get_genotype_matrix = function(genofile, snp_ids, sample_ids) {
    message("Loading complete genotype matrix for cis SNPs...")
    
    # Get genotypes for all cis SNPs at once
    geno_mat = snpgdsGetGeno(
        genofile,
        sample.id = sample_ids,
        snp.id = snp_ids,
        with.id = TRUE,
        verbose = FALSE
    )
    
    return(as.data.frame(geno_mat$genotype))
}

#' Process chromosome data for a specific brain region
#' @param chr Chromosome number
#' @param geno_data Full genotype matrix
#' @param snp_info SNP information
#' @param config Configuration
#' @param region Brain region being processed
process_chromosome = function(chr, geno_data, snp_info, config, region) {
    message(sprintf("Processing chromosome %d for %s", chr, region))
    
    # Filter SNPs for chromosome
    chr_snps = snp_info$SNP[snp_info$Chromosome == chr]
    
    if(length(chr_snps) == 0) {
        message(sprintf("No SNPs found for chromosome %d in %s", chr, region))
        return(NULL)
    }
    
    ## Subset genotype 
    genotype = as.data.frame(geno_data[, which(colnames(geno_data) %in% chr_snps)])
    ## Remove SNPs with missing values (Elastic Net training does not handle NAs)
    genotype = genotype[, apply(genotype, 2, function(x) all(!is.na(x)))]
    
    # Create chromosome-specific data
    chr_data = list(
        genotype = genotype,
        snp_info = snp_info[which(snp_info$SNP %in% colnames(genotype)), ]
    )
    
    # Create region-specific output directory
    region_output_dir = file.path(config$output_dir, region)
    dir.create(region_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save to RDS
    output_file = file.path(region_output_dir, sprintf("chr%d_genotype.rds", chr))
    saveRDS(chr_data, file = output_file)
    
    message(sprintf("Saved %s", output_file))
}

#################################################
## Main Processing
#################################################

main = function(config) {
    # Create output directory
    dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Process each brain region
    for (region in config$brain_regions) {
        message(sprintf("\nProcessing brain region: %s", region))
        
        # Create region-specific directory
        region_dir = file.path(config$output_dir, region)
        dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)
        
        message("Loading gene positions...")
        gene_positions = get_gene_positions(config$rnaseq, region)
        
        message("Processing genotype data")
        
        # Convert PLINK to GDS if needed
        gds_file = paste0(config$genotype_path, ".gds")
        if(!file.exists(gds_file)) {
            snpgdsBED2GDS(
                paste0(config$genotype_path, ".bed"),
                paste0(config$genotype_path, ".fam"),
                paste0(config$genotype_path, ".bim"),
                paste0(config$genotype_path, ".gds")
            )
        }
        
        # Open GDS file
        genofile = snpgdsOpen(gds_file, readonly = TRUE, allow.fork = TRUE)
        on.exit({
          snpgdsClose(genofile)
          file.remove(gds_file)
            
        })
        
        # Get SNP information and find cis SNPs
        message("Finding SNPs in cis windows...")
        snp_info = snpgdsSNPList(genofile)
        cis_snps = find_cis_snps(snp_info, gene_positions, config$cis_window)
        
        message(sprintf("Found %d SNPs in cis windows for %s", length(cis_snps), region))
        
        # Filter SNP info for cis SNPs only
        snp_info = snp_info[snp_info$snp.id %in% cis_snps, ]
        
        # Create SNP map
        snp_map = data.frame(
            SNP        = as.character(snp_info$snp.id),
            Chromosome = snp_info$chromosome,
            Position   = snp_info$position,
            Allele     = as.character(snp_info$allele),
            stringsAsFactors = FALSE
        )
        
        # Process allele information
        s                 = strsplit(snp_map$Allele, split = "/")
        snp_map$Al1       = sapply(s, `[`, 1)
        snp_map$Al2       = sapply(s, `[`, 2)
        snp_map$chr_pos   = paste0("chr", snp_map$Chromosome, ":", snp_map$Position)
        snp_map$varID     = paste0(snp_map$chr_pos, ":", snp_map$Al1, ":", snp_map$Al2)
        rownames(snp_map) = snp_map$SNP
        
        # Save SNP map for this region
        snp_map_file = file.path(region_dir, "snp_map.rds")
        saveRDS(snp_map, file = snp_map_file)
        message(sprintf("Saved SNP map to %s", snp_map_file))
        
        # Get sample IDs
        rna_seq    = get(load(config$rnaseq))
        sample_ids = rownames(rna_seq[[region]]$colData) 
        
        # Get complete genotype matrix for cis SNPs
        geno_data = get_genotype_matrix(
            genofile   = genofile,
            snp_ids    = cis_snps,
            sample_ids = sample_ids
        )
        
        # Set column names for the genotype matrix
        rownames(geno_data) = sample_ids
        colnames(geno_data) = cis_snps
        
        # Process each chromosome for this region
        chromosomes = sort(unique(snp_map$Chromosome))
        for(chr in chromosomes) {
            process_chromosome(chr, geno_data = geno_data, snp_info = snp_map, config, region)
        }
        
        message(sprintf("Completed processing for %s", region))
    }
    
    message("All brain regions processing complete!")
}

# Run the pipeline
tryCatch({
    main(CONFIG)
}, error = function(e) {
    message("Error in processing: ", e$message)
}, finally = {
    gc()
}) 
