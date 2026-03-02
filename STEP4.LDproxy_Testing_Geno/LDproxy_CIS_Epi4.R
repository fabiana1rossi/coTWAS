library(data.table)
library(parallel)
library(dplyr)
library(purrr)
library(limma)

# Set working directory to project root
# setwd("path/to/project/root")


#####
### Params
####

model="CIS_EpiXcan"
model_cis = "CIS"
model_epi = "EpiXcan"

### Params
target_file_path = c("./genotype/")
site               = "your_genotype.bim" # GTeX.v9.phased.geno.maf.bfile.updated.EUR    CMC_EUR.IBD.geno0.02.hwe10e6.mind0.02.maf0.01.hg20  LIBD_TopMed.geno.maf.EUR.chr1_22.noindels.overlap.GTEX
name_site          = "your_genotype_name"     # GTeX.v9.phased.geno.maf.bfile.updated.EUR   CMC_EUR.IBD.geno0.02.hwe10e6.mind0.02.maf0.01.hg20  LIBD_TopMed.geno.maf.EUR.chr1_22.noindels.overlap.GTEX


## Extract target genotype SNP bim files
bim_file        = fread(file.path(target_file_path,site))


## Get the LD table for the target dataset
load(paste0("./STEP4.LDproxy_Testing_Geno/LD_output/",model,".", name_site, "_LDtable.RData")) ## from previous script



########################
####### Add SNP annotation
#####################
print("Getting SNPs info")

all_weights_cis = readRDS("./STEP2.Training_CIS/output/database/weights_LD.Rds") ## dataframe containing gene-snps-ref_allele-alt_allele-POS-CHR  (from SQLite DBs)
all_weights_epi = readRDS("./STEP2.Training_CIS/output_EpiXcan/database/weights_LD.Rds") ## dataframe containing gene-snps-ref_allele-alt_allele-POS-CHR  (from SQLite DBs)

## Extract SNP map files
SNPmap_df      = bim_file


## Check alleles
invert_alleles = function(ref,target.LD,ref_allele.model, eff_allele.model,check,chr,POS,UNPHASED_R,MAJ_A,MAJ_B,REF.g1000,ALT.g1000){
  
  if(eff_allele.model == MAJ_A & sign(UNPHASED_R)>0){
    if(REF.g1000==MAJ_B){ return(list(ref=ref,target.LD=target.LD,ref_allele.model=ALT.g1000, eff_allele.model= MAJ_B,check=check,chr=as.character(chr),POS=as.character(POS),UNPHASED_R=as.numeric(UNPHASED_R),MAJ_A=MAJ_A,MAJ_B=MAJ_B,REF.g1000=REF.g1000,ALT.g1000=ALT.g1000))}
    if(REF.g1000!=MAJ_B){ return(list(ref=ref,target.LD=target.LD,ref_allele.model=REF.g1000, eff_allele.model= MAJ_B,check=check,chr=as.character(chr),POS=as.character(POS),UNPHASED_R=as.numeric(UNPHASED_R),MAJ_A=MAJ_A,MAJ_B=MAJ_B,REF.g1000=REF.g1000,ALT.g1000=ALT.g1000))}
  }
  if(eff_allele.model == MAJ_A & sign(UNPHASED_R)<0){
    if(REF.g1000==MAJ_B){ return(list(ref=ref,target.LD=target.LD,ref_allele.model=MAJ_B, eff_allele.model=ALT.g1000,check=check,chr=as.character(chr),POS=as.character(POS),UNPHASED_R=as.numeric(UNPHASED_R),MAJ_A=MAJ_A,MAJ_B=MAJ_B,REF.g1000=REF.g1000,ALT.g1000=ALT.g1000 ))}
    if(REF.g1000!=MAJ_B){ return(list(ref=ref,target.LD=target.LD,ref_allele.model=MAJ_B, eff_allele.model= REF.g1000,check=check,chr=as.character(chr),POS=as.character(POS),UNPHASED_R=as.numeric(UNPHASED_R),MAJ_A=MAJ_A,MAJ_B=MAJ_B,REF.g1000=REF.g1000,ALT.g1000=ALT.g1000))}
    
  }
  if(eff_allele.model != MAJ_A & sign(UNPHASED_R)>0){
    if(REF.g1000==MAJ_B){ return(list(ref=ref,target.LD=target.LD,ref_allele.model=MAJ_B, eff_allele.model= ALT.g1000,check=check,chr=as.character(chr),POS=as.character(POS),UNPHASED_R=as.numeric(UNPHASED_R),MAJ_A=MAJ_A,MAJ_B=MAJ_B,REF.g1000=REF.g1000,ALT.g1000=ALT.g1000))}
    if(REF.g1000!=MAJ_B){ return(list(ref=ref,target.LD=target.LD,ref_allele.model=MAJ_B, eff_allele.model= REF.g1000,check=check,chr=as.character(chr),POS=as.character(POS),UNPHASED_R=as.numeric(UNPHASED_R),MAJ_A=MAJ_A,MAJ_B=MAJ_B,REF.g1000=REF.g1000,ALT.g1000=ALT.g1000))}
  }
  if(eff_allele.model != MAJ_A & sign(UNPHASED_R)<0){
    if(REF.g1000==MAJ_B){ return(list(ref=ref,target.LD=target.LD,ref_allele.model=ALT.g1000, eff_allele.model= MAJ_B,check=check,chr=as.character(chr),POS=as.character(POS),UNPHASED_R=as.numeric(UNPHASED_R),MAJ_A=MAJ_A,MAJ_B=MAJ_B,REF.g1000=REF.g1000,ALT.g1000=ALT.g1000))}
    if(REF.g1000!=MAJ_B){ return(list(ref=ref,target.LD=target.LD,ref_allele.model=REF.g1000, eff_allele.model= ALT.g1000,check=check,chr=as.character(chr),POS=as.character(POS),UNPHASED_R=as.numeric(UNPHASED_R),MAJ_A=MAJ_A,MAJ_B=MAJ_B,REF.g1000=REF.g1000,ALT.g1000=ALT.g1000 ))}
  }
}


## subset with specific dataset
LD.SNPs.bim = inner_join(LD, SNPmap_df, by=c("target.LD"="V2")) %>% dplyr::select(-c(V1,V3,V5,V6)) %>% dplyr::rename(POS = V4)

## CIS
LD.SNPs.bim_cis = left_join(LD.SNPs.bim, all_weights_cis, by=c("ref"),relationship = "many-to-many") %>% dplyr::rename(chr = CHR)
colnames(LD.SNPs.bim_cis)[grepl("REF$|EFF$", colnames(LD.SNPs.bim_cis))] =c("ref_allele.model","eff_allele.model")
LD.SNPs.bim_cis = LD.SNPs.bim_cis %>% dplyr::select(ref,target.LD,ref_allele.model, eff_allele.model,check,chr,POS,UNPHASED_R,MAJ_A,MAJ_B, REF.g1000,ALT.g1000)
LD.SNPs.bim_cis = na.omit(LD.SNPs.bim_cis)

print(paste0("inverting alleles for site: ",name_site))
index=which(LD.SNPs.bim_cis$check=="different")
LD.SNPs.bim_cis$chr = as.character(LD.SNPs.bim_cis$chr)
LD.SNPs.bim_cis$POS = as.character(LD.SNPs.bim_cis$POS)
LD.SNPs.bim_cis[index,] = pmap_dfr(LD.SNPs.bim_cis[index,], invert_alleles)


LD.SNPs.bim_cis = LD.SNPs.bim_cis %>% distinct() %>% mutate(Site=name_site) %>%
  dplyr::select(ref,target.LD,check,chr,POS, ref_allele.model,eff_allele.model,Site) %>%
  dplyr::rename(chrom_start="POS",ref_allele=ref_allele.model,eff_allele=eff_allele.model)

## EpiXcan
LD.SNPs.bim_epi = left_join(LD.SNPs.bim, all_weights_epi, by=c("ref"),relationship = "many-to-many") %>% dplyr::rename(chr = CHR)
colnames(LD.SNPs.bim_epi)[grepl("REF$|EFF$", colnames(LD.SNPs.bim_epi))] =c("ref_allele.model","eff_allele.model")
LD.SNPs.bim_epi = LD.SNPs.bim_epi %>% dplyr::select(ref,target.LD,ref_allele.model, eff_allele.model,check,chr,POS,UNPHASED_R,MAJ_A,MAJ_B, REF.g1000,ALT.g1000)
LD.SNPs.bim_epi = na.omit(LD.SNPs.bim_epi)

print(paste0("inverting alleles for site: ",name_site))
index=which(LD.SNPs.bim_epi$check=="different")
LD.SNPs.bim_epi$chr = as.character(LD.SNPs.bim_epi$chr)
LD.SNPs.bim_epi$POS = as.character(LD.SNPs.bim_epi$POS)
LD.SNPs.bim_epi[index,] = pmap_dfr(LD.SNPs.bim_epi[index,], invert_alleles)


LD.SNPs.bim_epi = LD.SNPs.bim_epi %>% distinct() %>% mutate(Site=name_site) %>%
  dplyr::select(ref,target.LD,check,chr,POS, ref_allele.model,eff_allele.model,Site) %>%
  dplyr::rename(chrom_start="POS",ref_allele=ref_allele.model,eff_allele=eff_allele.model)


## Save
LD = distinct(LD.SNPs.bim_cis)
save(LD, file=paste0("./STEP4.LDproxy_Testing_Geno/LD_output/",model_cis,".", name_site, "_LDtable.RData"))

LD = distinct(LD.SNPs.bim_epi)
save(LD, file=paste0("./STEP4.LDproxy_Testing_Geno/LD_output/",model_epi,".", name_site, "_LDtable.RData"))

