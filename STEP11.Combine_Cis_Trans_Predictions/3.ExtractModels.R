# Set working directory to project root
# setwd("path/to/project/root")


tissues = c("dlpfc")

full = sapply(tissues,function(x){
  
  if(x!="amygdala"){
    epixcan = get(load(paste0("./predictions/EpiXcan/EpiXcan_",x,"_testing_dataset_name_performance_selected.RData")))
    epixcan = epixcan$gene
    
    module=get(load(paste0("./predictions/MODULE/MODULE_averaged/MODULE_replicable_testing_dataset_name_performance_",x,".RData")))
    module = module$gene[module$network %in% "Averaged_Ensemble"]
    
    ingene= get(load(paste0("./predictions/INGENE/INGENE_averaged/INGENE_replicable_testing_dataset_name_performance_",x,".RData")))
    ingene = ingene$gene[ingene$network %in% "Averaged_Ensemble"]
    
    cis=get(load(paste0("./predictions/CIS/CIS_",x,"_testing_dataset_name_performance_selected.RData")))
    cis = cis$gene
    
    
    unique.epixcan = data.frame(gene = setdiff(epixcan,union(union(module, ingene),cis)), model = "E") 
    unique.module = data.frame(gene = setdiff(module,union(union(epixcan, ingene),cis)),  model = "M")
    unique.ingene = data.frame(gene = setdiff(ingene,union(union(module, epixcan),cis)),  model = "I")
    unique.cis    = data.frame(gene = setdiff(cis,union(union(module, epixcan),ingene)),  model    = "C")
    unique = rbind(unique.epixcan,unique.module,unique.ingene,unique.cis)
    
    common.epixcan = epixcan[!epixcan %in% unique.epixcan$gene]
    common.module  = module[!module %in% unique.module$gene]
    common.ingene  = ingene[!ingene %in% unique.ingene$gene]
    common.cis     = cis[!cis %in% unique.cis$gene]
    
    ## Combination of 3-4 models
    common.EMI = data.frame(gene = Reduce(intersect,list(common.epixcan,common.ingene,common.module)), model = "EMI")
    common_EMI = common.EMI$gene
    
    common.MIC = data.frame(gene = Reduce(intersect,list(common.cis,common.ingene,common.module)), model = "MIC")
    common_MIC = common.MIC$gene
    
    ## Remaining genes
    
    ## Containing EpiXcan
    remaining_genes_epi = setdiff(common.epixcan, common_EMI)
    
    # Combination of 2 models containing epixcan
    common.EI.EM = data.frame(
      gene = remaining_genes_epi,
      model = ifelse(
        remaining_genes_epi %in% module & !remaining_genes_epi  %in% ingene, "EM", "EI")
    )
    common_EI.EM = common.EI.EM$gene
    
    ## Containing MODULE
    remaining_genes_module = setdiff(common.module, Reduce(union,list(common_EMI,common_MIC)))
    # Combination of 2 models containing module
    common.MI.MC = data.frame(
      gene = remaining_genes_module,
      model = ifelse(
        remaining_genes_module %in% ingene & !remaining_genes_module  %in% c(cis,epixcan), "MI",
        ifelse(remaining_genes_module %in% cis & !remaining_genes_module  %in% c(ingene,epixcan), "MC","DUPLICATED"))
    )
    common.MI.MC = common.MI.MC %>% filter(!model=="DUPLICATED")
    common_MI.MC = common.MI.MC$gene
    
    ## Containing INGENE
    remaining_genes_ingene = setdiff(common.ingene, Reduce(union,list(common_EMI,common_MIC)))
    # Combination of 2 models containing ingene
    common.IC = data.frame(
      gene = remaining_genes_ingene,
      model = ifelse(
        remaining_genes_ingene %in% cis & !remaining_genes_ingene  %in% c(module,epixcan), "IC", "DUPLICATED")
    )
    
    common.IC = common.IC %>% filter(!model=="DUPLICATED")
    common_IC  = common.IC$gene
    
    
    
    ## bind
    common = rbind(common.EMI,common.MIC,common.EI.EM,common.MI.MC,common.IC)
    
  }else{
    
  module=get(load(paste0("./predictions/MODULE/MODULE_averaged/MODULE_replicable_testing_dataset_name_performance_",x,".RData")))
  module = module$gene[module$network %in% "Averaged_Ensemble"]
  
  ingene= get(load(paste0("./predictions/INGENE/INGENE_averaged/INGENE_replicable_testing_dataset_name_performance_",x,".RData")))
  ingene = ingene$gene[ingene$network %in% "Averaged_Ensemble"]
  
  cis=get(load(paste0("./predictions/CIS/CIS_",x,"_testing_dataset_name_performance_selected.RData")))
  cis = cis$gene
  
  
  unique.module = data.frame(gene = setdiff(module,union(ingene,cis)),  model = "M")
  unique.ingene = data.frame(gene = setdiff(ingene,union(module,cis)),  model = "I")
  unique.cis    = data.frame(gene = setdiff(cis,union(module,ingene)),  model    = "C")
  unique = rbind(unique.module,unique.ingene,unique.cis)
  
  
  common.module  = module[!module %in% unique.module$gene]
  common.ingene  = ingene[!ingene %in% unique.ingene$gene]
  common.cis     = cis[!cis %in% unique.cis$gene]
  
  ## Combination of 3-4 models
  common.MIC = data.frame(gene = Reduce(intersect,list(common.cis,common.ingene,common.module)), model = "MIC")
  common_MIC = common.MIC$gene
  
  ## Remaining genes
  
  ## Containing MODULE
  remaining_genes_module = setdiff(common.module, Reduce(union,list(common_MIC)))
  # Combination of 2 models containing module
  common.MI.MC = data.frame(
    gene = remaining_genes_module,
    model = ifelse(
      remaining_genes_module %in% ingene & !remaining_genes_module  %in% c(cis), "MI",
      ifelse(remaining_genes_module %in% cis & !remaining_genes_module  %in% c(ingene), "MC","DUPLICATED"))
  )
  common.MI.MC = common.MI.MC %>% filter(!model=="DUPLICATED")
  common_MI.MC = common.MI.MC$gene
  
  ## Containing INGENE
  remaining_genes_ingene = setdiff(common.ingene, Reduce(union,list(common_MIC)))
  # Combination of 2 models containing ingene
  common.IC = data.frame(
    gene = remaining_genes_ingene,
    model = ifelse(
      remaining_genes_ingene %in% cis & !remaining_genes_ingene  %in% c(module), "IC", "DUPLICATED")
  )
  
  common.IC = common.IC %>% filter(!model=="DUPLICATED")
  common_IC  = common.IC$gene
  
  
  
  ## bind
  common = rbind(common.MIC,common.MI.MC,common.IC)
  }
  
  
  dplyr::lst(unique,common)
},USE.NAMES = T,simplify = F)

full = purrr::transpose(full)
save(full,file = "./predictions/Combined/EMIC_Common_Unique_updated.RData")

Common_Genes = full[2]
Common_Genes$Chosen.genes = lapply(Common_Genes$common,function(x) x$gene)
save(Common_Genes,file = "./predictions/Combined/EMIC_Common.Genes.Models_updated.RData")

Unique_Genes = full[1]
Unique_Genes$Chosen.genes = lapply(Unique_Genes$unique,function(x) x$gene)
save(Unique_Genes,file = "./predictions/Combined/EMIC_Unique.Genes.Models_updated.RData")

