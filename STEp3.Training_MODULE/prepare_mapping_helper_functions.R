

######################
## 1. Process SNPs ##
######################

process_snps = function(module_file, snplist_file, expr_file, geno_dir,geno_name, output_dir, brain_regions, n_cores) {
  
 ## Load data
  networkList    = readRDS(module_file)     # Load coexpression networks
  snpList        = get(load(snplist_file)) # loads snpList.100kbp
  training_expr  = get(load(expr_file))    # loads RNAseq training
  ## Subset expression data object with avaiable brain regions
  training_expr  = training_expr[brain_regions]
  
  # Get gene sets for each brain region
  gene_sets = map(training_expr, ~{
    colnames(.x[["assays"]][["expression"]])
  })
  
  # Get  SNPs annotated to  brain region genes
  snpList_brain_region = map(names(gene_sets) %>% set_names(.,.), ~{
    genes = gene_sets[[.x]]
    snpList[which(names(snpList) %in% genes)]
  })
  
  # Subset SNPs annotated to brain regions coexpression modules also present in genotype (bim) 
  snpList_modules = map(names(snpList_brain_region), ~{
    brain_region = .x
    
    # Get genes in network/modules
    network_genes   = unique(unlist(networkList))
    sublist_mod     = snpList_brain_region[[brain_region]][which(names(snpList_brain_region[[brain_region]]) %in% network_genes)]
    
    # Filter SNPs present in genotype data
    bim = read.delim(paste0(geno_dir, 
                            geno_name,'.bim'), header=F)
    
    sublist = mclapply(1:length(sublist_mod), function(i) {
      print(paste0(brain_region, ': ', i, '/',length(sublist_mod)))
      sublist_mod[[i]][which(sublist_mod[[i]] %in% bim$V2)]
    }, mc.cores = n_cores)
    
    names(sublist) = names(sublist_mod)
    save(sublist, file=paste0(output_dir, "snpList_",brain_region,".RData"))
    
    return(sublist)
  })
  
  return(snpList_modules)
}

