#' Load required packages
load_dependencies = function() {
    # Install and load CRAN packages
    cran_packages = c("purrr", "dplyr", "parallel")
    for(pkg in cran_packages) {
        if(!require(pkg, character.only = TRUE)) {
            install.packages(pkg)
            library(pkg, character.only = TRUE)
        }
    }

    # Install and load Bioconductor packages
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    bioc_packages = c("limma", "SummarizedExperiment") 
    for(pkg in bioc_packages) {
        if(!require(pkg, character.only = TRUE)) {
            BiocManager::install(pkg)
            library(pkg, character.only = TRUE)
        }
    }
}


#' Filter module list based on configuration
filter_networks = function(networks, config) {
    networks   = networks[!grepl(config$excluded_networks,names(networks))]
    networks   = map(names(networks) %>% set_names(.,.), function(modules){
      modules   = networks[[modules]][!grepl(config$excluded_modules, names( networks[[modules]]))]
     return(modules)
  })
    
   return(networks)
}



#' Process modules in parallel 
#' @param module_data List containing all modules in network and data
#' @return List of processed module data for each module
process_module = function(module_data, n_cores,brain_region, network, output_path) {
    
    
    # Process each module in parallel using mclapply
    results = mclapply(names(module_data$modules), function(module) {
        tryCatch({
            # Get SNP list for the brain region from module_data
             # Filter SNPMap based on config$snp_list specific for brain region
            snpList_region = module_data$snps 
            ## Subset snp annotated to genes in region, with genes in module
            snpList_module = snpList_region[names(snpList_region) %in% module_data$modules[[module]]]
          
            # Extract module expression data and subset with module genes
            mod_expr = module_data$ref_matrix[, which(colnames(module_data$ref_matrix) %in% module_data$modules[[module]]), drop=FALSE]
            
            ## Match samples with genotype
            mod_expr         = mod_expr[which(rownames(mod_expr) %in% rownames(module_data$geno)),]
            module_data$geno = module_data$geno[which(rownames(module_data$geno) %in% rownames(mod_expr)),]
            
            mod_expr = mod_expr[match(rownames(module_data$geno),rownames(mod_expr)),]
            
            ## if only 1 gene
            if(ncol(mod_expr) <= 1) return(NULL)
            
            # Calculate PC1 
            pc1 = prcomp(mod_expr, scale=TRUE, center=TRUE)$x[,1]
            
            # Subset genotype data with SNPs annotated to module 
            geno_mod = module_data$geno[, unlist(unique(snpList_module)), drop=FALSE]
            map_mod  = module_data$map[unlist(unique(snpList_module)),]
            
            # Return results
            module_res = list(
                eigen0    = pc1,
                genes     = module_data$modules[[module]],
                snps      = unlist(unique(snpList_module)),
                geno.mod  = geno_mod,
                map.mod   = map_mod
            )
            
            ## Save module list
            save(module_res, file = sprintf("%s%s_network_%s_module_%s.RData", output_path, brain_region, network, module))
        }, error = function(e) {
            message("Error processing module: ", e$message)
            return(NULL)
        })
    }, mc.cores = n_cores)
    
    # Name results with module names
    names(results) = names(module_data$modules)
    
    return(results)
}


#' Load and process genotype data
#' @param config Configuration list containing paths and parameters
#' @param brain_region String specifying the brain region
#' @param expr Expression matrix with samples as rows
#' @return List containing filtered genotype matrix and SNP map
load_genotype_data = function(config, brain_region, expr_samples) {
    # genotype base filename
    base_file =  config$genotype_path
 
    ## open FAM to match expression samples
    fam          = read.delim(paste0(base_file, ".fam"),header = F)
    expr_samples = expr_samples[which(expr_samples %in% fam$V2)]
    
    # Convert PLINK to GDS format
    snpgdsBED2GDS(
        paste0(base_file, ".bed"),
        paste0(base_file, ".fam"),
        paste0(base_file, ".bim"),
        paste0(base_file, ".gds")
    )
    
    # Open GDS file
    genofile = snpgdsOpen(paste0(base_file, ".gds"),
                         readonly = FALSE,
                         allow.duplicate = FALSE,
                         allow.fork = TRUE)
    
    # Create SNP map
    options(stringsAsFactors = FALSE)
    snpMap = snpgdsSNPList(genofile, sample.id = expr_samples)
    snpMap.df = data.frame(
        SNP = as.character(snpMap$snp.id),
        Chromosome = snpMap$chromosome,
        Position = snpMap$position,
        Allele = as.character(snpMap$allele),
        Freq = snpMap$afreq
    )
    
    # Process allele information
    s = strsplit(snpMap.df$Allele, split = "/")
    snpMap.df$Al1 = sapply(s, function(i) i[1])
    snpMap.df$Al2 = sapply(s, function(i) i[2])
    snpMap.df$chr_pos = paste("chr", snpMap.df$Chromosome, snpMap.df$Position, sep = ":")
    snpMap.df$varID = paste0("chr", snpMap.df$Chromosome, ":", snpMap.df$Position, 
                            ":", snpMap.df$Al1, ":", snpMap.df$Al2)
    rownames(snpMap.df) = snpMap.df$SNP
    
    # Filter SNPMap based on config$snp_list specific for brain region
    snpList_region = get(load(config$snplist_file[grepl(brain_region,config$snplist_file)])) 
    
    valid_snps = intersect(snpMap.df$SNP, unique(unlist(snpList_region)))
    snpMap.df = snpMap.df[valid_snps, ]
    
    # Get genotype data for filtered SNPs
    g = snpgdsGetGeno(genofile,
                      sample.id = expr_samples,
                      snp.id = valid_snps,
                      with.id = TRUE,
                      verbose = TRUE)
    
    # Process genotype matrix
    geno = as.data.frame(g$genotype)
    colnames(geno) = g$snp.id
    rownames(geno) = g$sample.id
    
    
    # Clean up
    snpgdsClose(genofile)
    file.remove(paste0(base_file, ".gds"))
    
    return(list(
        geno.test = geno,
        map.test = snpMap.df
    ))
}


#' Process a single brain region
#' @param brain_data Brain region specific data
#' @param config Pipeline configuration
#' @param brain_region Current brain region name
#' @return NULL (saves results to files)
process_brain_region = function(brain_data, config, brain_region) {
    message(sprintf("Processing brain region: %s", brain_region))
    
    ## Get SNPs annotated to genes in region 
    snpList_region = get(load(config$snplist_file[grepl(brain_region,config$snplist_file)])) 
  
    # Get genotype data
    expr_samples = brain_data$colData$BrNum  ## get brain specific samples from expression
    geno_data    = load_genotype_data(config, brain_region, expr_samples)
    
   
    # Process each network and module in parallel
    mclapply(names(config$networks), function(network_name) {
        message(sprintf("Processing network: %s", network_name))
        
        module_data = list(
            ref_matrix = as.data.frame(brain_data$assays$expression),
            modules    = config$networks[[network_name]],
            snps       = snpList_region,
            geno       = geno_data[[1]],
            map        = geno_data[[2]]
        )
        
        ## process and save a file for each module
        process_module(module_data, n_cores = config$mc_cores, brain_region, network_name, config$output_dir) ## Specify more cores to run code in parallel on ubuntu terminal
        
        # if(!is.null(res_network)) {
        #     output_file = file.path(
        #         config$output_dir,
        #         paste0(brain_region, '_network_', network_name,
        #               '.RData')
        #     )
            #save(res_network, file = output_file)
        #}
    }, mc.cores = config$mc_cores)
}


############################
#' Main pipeline function
#######################
#' Parameters:
#' config: List containing the following fields:
#'   - expression_file: Path to RData file with expression data
#'   - module_list_file: Path to RDS file with gene modules
#'   - genotype_path: Directory containing genotype files
#'   - output_dir: Directory for output files
#'   - brain_regions: Vector of brain region names 
#'   - snp_list: List of SNPs to analyze
#'   - ethnicity: String specifying ethnicity for file naming
#'   - excluded_patterns: Regex pattern for modules to exclude
#'   - excluded_modules: Vector of module names to exclude
#'   - mc_cores: Number of cores for parallel processing

run_mapping_pipeline = function(config) {
  # Load dependencies
  load_dependencies()
  
  # Load and filter modules
  networks        = readRDS(config$networks_list_file)
  config$networks = filter_networks(networks, config)
  
  # Load expression data
  rnaseq  = get(load(config$expression_file))
  
  # Process specified brain regions
  brain_regions = config$brain_regions
  
  # Process each brain region
  for(brain_region in brain_regions) {
    process_brain_region(
      rnaseq[[brain_region]],
      config,
      brain_region
    )
  }
}

