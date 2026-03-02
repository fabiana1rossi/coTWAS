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

### Params
target_file_path = c("./genotype/")
site               = "your_geno_bim"
name_site          = "your_geno_name"


## Extract target genotype SNP bim files
bim_file        = fread(file.path(target_file_path,site))



## Get the LD table for the target dataset
load(paste0("./STEP4.LDproxy_Testing_Geno/LD_output/MODULE.", name_site, "_LDtable.RData"))



########################
####### Add SNP annotation
#####################
print("Getting SNPs info")

all_weights = get(load("./STEP3.Training_MODULE/output/database/all_weights_LD.RData"))
colnames(all_weights)[c(1,3)] = c("ref","pos")

## Extract SNP map files
SNPmap_df      = bim_file

## subset with specific dataset
LD.SNPs.bim = inner_join(LD, SNPmap_df, by=c("target.LD"="V2")) %>% dplyr::select(-c(V3,V5,V6)) %>% dplyr::rename(chr = V1, POS = V4)
LD.SNPs.bim = left_join(LD.SNPs.bim, all_weights, by=c("ref"),relationship = "many-to-many")
colnames(LD.SNPs.bim)[grepl("REF$|EFF$", colnames(LD.SNPs.bim))] =c("ref_allele.model","eff_allele.model")

LD.SNPs.bim = LD.SNPs.bim %>% dplyr::select(ref,target.LD,ref_allele.model, eff_allele.model,check,chr,POS,UNPHASED_R,MAJ_A,MAJ_B, REF.g1000,ALT.g1000)

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


print(paste0("inverting alleles for site: ",name_site))
index=which(LD.SNPs.bim$check=="different")
LD.SNPs.bim$chr = as.character(LD.SNPs.bim$chr)
LD.SNPs.bim$POS = as.character(LD.SNPs.bim$POS)
LD.SNPs.bim[index,] = pmap_dfr(LD.SNPs.bim[index,], invert_alleles)


LD.SNPs.bim = LD.SNPs.bim %>% distinct() %>% mutate(Site=name_site) %>%
dplyr::select(ref,target.LD,check,chr,POS, ref_allele.model,eff_allele.model,Site) %>%
dplyr::rename(chrom_start="POS",ref_allele=ref_allele.model,eff_allele=eff_allele.model)



LD = distinct(LD.SNPs.bim)


save(LD, file=paste0("./STEP4.LDproxy_Testing_Geno/LD_output/MODULE.", name_site, "_LDtable.RData"))


