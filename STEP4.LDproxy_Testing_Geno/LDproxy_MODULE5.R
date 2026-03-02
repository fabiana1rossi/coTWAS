library(RSQLite)
library(purrr)
library(dplyr)
library(limma)
library(parallel)
library(DBI)


# Set working directory to project root
# setwd("path/to/project/root")


###############
## PARAMS
###############
models          = c('MODULE')
target_file_path = c("./genotype/")
name_site          = "your_geno_name" #CMC_EUR.IBD.geno0.02.hwe10e6.mind0.02.maf0.01.hg20




#################
## FUNCTION
#################
# Function to process each gene in the PRedictDB
process_weights = function(weights_df, LD.subset_df, model,gene_annot) {
  
  
  matching_snps        = weights_df$snp[weights_df$snp %in% LD.subset_df$ref]
  weights_subset       = weights_df[weights_df$snp %in% matching_snps,]
  colnames(weights_subset)[1:6] = c("gene","rsid","varID","ref_allele","eff_allele","weight")
  weights_nomatching   = weights_df[-which(weights_df$snp %in% matching_snps),]
  weights_nomatching$Site = weights_nomatching$original.rsid = NA
  
  if(any(grepl("network",colnames(weights_nomatching)))){
    weights_nomatching  = weights_nomatching %>% dplyr::select(gene,network,snp,varID,ref_vcf,alt_vcf,weights,Site,original.rsid) %>% dplyr::rename(rsid = "snp", ref_allele = "ref_vcf", eff_allele = "alt_vcf", weight = "weights")
  }else{
    weights_nomatching  = weights_nomatching %>% dplyr::select(gene,snp,varID,ref_vcf,alt_vcf,weights,Site,original.rsid)%>% dplyr::rename(rsid = "snp", ref_allele = "ref_vcf", eff_allele = "alt_vcf", weight = "weights")
  }
  
  LD.subset_df         = LD.subset_df[LD.subset_df$ref %in% matching_snps,]
  
  
  ##take the most correlated one per ref --> si perdono delle snp ma la computazioen è piu veloce
  LD.subset_df_nodup   = LD.subset_df[-which(duplicated(LD.subset_df$ref)),]
  ## add varID
  LD.subset_df_nodup$varID  = paste0("chr",LD.subset_df_nodup$chr,":",LD.subset_df_nodup$chrom_start,":",LD.subset_df_nodup$ref_allele,":",LD.subset_df_nodup$eff_allele)
  
  ## join dataframes
  weights_LD = inner_join(weights_subset, LD.subset_df_nodup, by=c("rsid"="ref")) 
  
  
  # ## Add gene position
  # weights_LD = inner_join(weights_LD,gene_annot, by=c("gene"="ensemblID"))
  
  ## Substitute new snp
  if(any(grepl("network",colnames(weights_LD)))){
    weights_LD           = weights_LD %>% dplyr::select(gene,network,target.LD,varID.y,ref_allele.y,eff_allele.y,weight,Site,rsid) 
    colnames(weights_LD) = c("gene","network","rsid","varID","ref_allele","eff_allele","weight","Site","original.rsid")
  }else{
    weights_LD           = weights_LD %>% dplyr::select(gene,target.LD,varID.y,ref_allele.y,eff_allele.y,weight,Site,rsid) 
    colnames(weights_LD) = c("gene","rsid","varID","ref_allele","eff_allele","weight","Site","original.rsid")
    
  }
  
  
  weights              = bind_rows(weights_LD, weights_nomatching) %>% distinct()

  return(weights)
}


#########################
## MAIN
#######################
start = Sys.time()
walk(models, function(model){
  regions    = c('dlpfc','acc','caudate','hippo','amygdala')
 
  target    = name_site
  
  ## get LD table
  LD = get(load(paste0("./STEP4.LDproxy_Testing_Geno/LD_output/",model,".", name_site, "_LDtable.RData")))
  
  ## Remove rows corresponding to LD SNPs with no annotation
  LD        = LD %>%
    filter(!(check == "different" & is.na(ref_allele)))
  gc()
  
 
  walk(regions, function(brain.region){
    print(paste0(model, " in ", target, '; ',toupper(brain.region)))
    
    ## get weights
    weights_region = readRDS(paste0("./STEP3.Training_MODULE/output/database/all_networks_weights_",brain.region,".rds"))
    
  
    if(file.exists(paste0("./PredictDB_LD/",model, "/" ,model,".",brain.region,".extras_weights_",name_site,".RData"))){return(NA)}      
      
    ## Match LD SNPs
    LD.subset     = LD %>% filter(check=='different',Site==name_site)
      
    gc()
      
      
    ## Subset of genes in bathces
    unique_snps =unique(weights_region$snp)
      
      
    batch_weights = process_weights(weights_df = weights_region, LD.subset_df = LD.subset, model = model, gene_annot = NULL)
    weights       = distinct(batch_weights)
      
      
    save(weights, file=paste0("./PredictDB_LD/",model, "/" ,model,".",brain.region,".extras_weights_",name_site,".RData"))
      
      
    print(paste0(name_site,' DONE'))

  }) # Close walk(regions, function(brain.region)
}) # Close walk(models, function(model)

print(Sys.time()-start)
