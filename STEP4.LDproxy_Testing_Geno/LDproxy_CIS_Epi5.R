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
models          = c('CIS','EpiXcan')
target_file_path = c("./genotype/")
name_site          = "your_genotype_name"     # GTeX.v9.phased.geno.maf.bfile.updated.EUR   CMC_EUR.IBD.geno0.02.hwe10e6.mind0.02.maf0.01.hg20  LIBD_TopMed.geno.maf.EUR.chr1_22.noindels.overlap.GTEX


#################
## FUNCTION
#################
# Function to process each gene in the PRedictDB
process_weights = function(weights_df, LD.subset_df, model,gene_annot) {
  
  matching_snps        = weights_df$rsid[weights_df$rsid %in% LD.subset_df$ref]
  weights_subset       = weights_df[weights_df$rsid %in% matching_snps,] %>% dplyr::select(gene, rsid,chr,pos, ref_allele, eff_allele,weight)
  colnames(weights_subset)= c("gene","rsid","chr","pos","ref_allele","eff_allele","weight")
  weights_nomatching   = weights_df[-which(weights_df$rsid %in% matching_snps),]
  weights_nomatching$Site = weights_nomatching$original.rsid = NA
  weights_nomatching  = weights_nomatching %>% dplyr::select(gene,rsid,chr,pos,ref_allele,eff_allele,weight,Site,original.rsid)
  
  
  LD.subset_df         = LD.subset_df[LD.subset_df$ref %in% matching_snps,]
  
  ##take the most correlated one per ref --> si perdono delle snp ma la computazioen è piu veloce
  LD.subset_df_nodup   = LD.subset_df[-which(duplicated(LD.subset_df$ref)),]

  ## join dataframes
  weights_LD = inner_join(weights_subset, LD.subset_df_nodup, by=c("rsid"="ref")) 
  
  ## Substitute new snp

  weights_LD           = weights_LD %>% dplyr::select(gene,target.LD,chr.x,pos,ref_allele.y,eff_allele.y,weight,Site,rsid) 
  colnames(weights_LD) = c("gene","rsid","chr","pos","ref_allele","eff_allele","weight","Site","original.rsid")
  
  
  
  
  weights              = bind_rows(weights_LD, weights_nomatching) %>% distinct()
  
  return(weights)
}


#########################
## MAIN
#######################
start = Sys.time()
walk(models, function(model){
  if(model=="CIS"){regions=c("dlpfc","caudate","hippo","amygdala","sACC","dACC"); file_path="output_LIBD"}else{regions=c("dlpfc","caudate","hippo","sACC","dACC");file_path="output_LIBD_EpiXcan"}
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
    weights_region = read.csv(paste0("./STEP2.Training_CIS/",file_path,"/database/weights_",brain.region,".csv"))
  
    
    if(file.exists(paste0("./PredictDB_LD/",model, "/" ,model,".",brain.region,".extras_weights_",name_site,".RData"))){return(NA)}      
    
    ## Match LD SNPs
    LD.subset     = LD %>% filter(check=='different',Site==name_site)
    
    gc()
    
    
    ## Subset of genes in bathces
    unique_snps = unique(weights_region$rsid)
    
    
    batch_weights = process_weights(weights_df = weights_region, LD.subset_df = LD.subset, model = model, gene_annot = NULL)
    weights       = distinct(batch_weights)
    
    
    save(weights, file=paste0("./PredictDB_LD/",model, "/" ,model,".",brain.region,".extras_weights_",name_site,".RData"))
    
    
    print(paste0(name_site,' DONE'))
    
  }) # Close walk(regions, function(brain.region)
}) # Close walk(models, function(model)

print(Sys.time()-start)

