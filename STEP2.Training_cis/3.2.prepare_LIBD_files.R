
# Set working directory - adjust path as needed
# setwd("path/to/STEP2.Training_cis/")

library(purrr)
library(dplyr)
library(data.table)
library(limma)

## Get unique snp_annot all chromosomes and convert genotypes files in RDAta
walk(c("dlpfc"), function(region){
  eqtls_all <- fread(sprintf("./output_LIBD/%s/LIBD_%s_all_chr_eqtls.txt",region,region), fill = TRUE)
  eqtls_all  = na.omit(eqtls_all)
 
  walk(c(1:22), function(chr){
   print(sprintf("Region %s, Chr %d", region,chr))
    geno = fread(sprintf("./genotype_by_chr_LIBD/%s/training_genotype.chr%d.txt",region,chr),sep = "\t")
    geno = as.data.frame(geno)
    rownames(geno) = geno[[1]]
    geno[[1]]  = NULL
    
   
    snp_annot = fread(sprintf("./genotype_by_chr_LIBD/%s/training_snp_annot.chr%d.txt",region,chr),sep = "\t")
    snp_annot = as.data.frame(snp_annot)
    snp_annot[[1]]  = NULL
    snp_annot  = snp_annot %>% dplyr::select(rsid, chr, pos,alt,ref,varID)
    snp_annot = distinct(snp_annot)
    rownames(snp_annot) = snp_annot$rsid
    snp_annot  = snp_annot %>% dplyr::rename(SNP="rsid", Chromosome="chr", Position="pos",Al1="alt",Al2="ref")
    snp_annot$Variant = paste0("chr",snp_annot$Chromosome,".",snp_annot$Position)
    snp_annot = snp_annot[which(rownames(snp_annot) %in% colnames(geno)),]
    
    
    saveRDS(list(genotype=geno, snp_info=snp_annot), file = sprintf("./genotype_by_chr_LIBD/%s/chr%d_genotype.Rds",region,chr))
    
    final_all = inner_join(snp_annot, eqtls_all, by ="Variant")
   
    
    saveRDS(final_all, file = sprintf("./output_LIBD/%s/chr%d_signifcant_ciseQTL.Rds",region,chr))
      
    
   
   
    })
 
})









