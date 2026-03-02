###############################################################################################
# Script to clean confounder effects from expression data
# This script performs:
# 1. Loads preprocessed RSE data with outliers removed
# 2. Handles duplicate samples (keeps RiboZeroGold or highest RIN PolyA samples)
# 3. Removes confounders using linear regression (Dataset, Protocol, Age, Sex, etc.)
# 4. Calculates and adds expression PCs to covariate model
# 5. Applies rank normalization and z-score transformation
###############################################################################################


library(jaffelab)
set.seed(123)
library(SummarizedExperiment)
library(purrr)

setwd('')


rse_preprocessed = readRDS("outlierRemoved.rds")        


###############
get.cleaned.obj = function(rse, rse.name, protect,  assay.to.clean){
   print(rse.name)
   rse$RIN = sapply(rse$RIN, mean)
   ## Check duplicated.
   print(paste0('duplicated BrNum: ',sum(duplicated(rse$BrNum))))
   BrNum_vector <- colData(rse) %>%
      as.data.frame() %>%
      dplyr::group_by(BrNum) %>%
      # Keep RiboZeroGold rows and PolyA rows with highest RIN 
      dplyr::slice(ifelse(
         all(Protocol == "PolyA" | Protocol == "RiboZeroGold" | Protocol == "RiboZeroHMR"),
         which.max(RIN),
         which(Protocol %in% c("RiboZeroGold"))
      )) %>%
      dplyr::ungroup() %>%
      dplyr::pull(BrNum)
   
   
   # Get the indices of the BrNum values in the SummarizedExperiment object
   indices <- match(BrNum_vector, rse$BrNum)
   # Subset the SummarizedExperiment object based on the indices
   rse <- rse[, indices]
   
   
   if(rse.name=='caudate'){
      rse = subset(rse,colData(rse)$Dx != 'MDD',subset = TRUE)  ## caudate has only 1 subject MDD
   }
   print(paste0('Number of Samples: ',length(rse$BrNum)))
   
   covs = as.data.frame(colData(rse))
   covs$Dataset = as.character(covs$Dataset)
   covs$Dataset = as.factor(covs$Dataset)
   covs$Protocol = as.character(covs$Protocol)
   covs$Protocol = factor(covs$Protocol)
   covs$Dx = as.character(covs$Dx)
   covs$Dx = as.factor(covs$Dx)

   exp = assays(rse)[[assay.to.clean]]
 
   ## Make mod. 
   if(grepl("dlpfc",rse.name)){
    mod = model.matrix(~   (Dx) + Age +  (Sex) + (mitoRate) + (rRNA_rate) + (totalAssignedGene) + (RIN)+
                            (Dataset) + C1 + C2 + C3 + C4+ C5, 
                       data = covs) ## Dataset & protocol colinear in DLPFC
   }
   if(grepl("hippo",rse.name)){
     mod = model.matrix(~   (Dx) + Age +  (Sex) + (mitoRate) + (rRNA_rate) + (totalAssignedGene) + (RIN)+(Dataset)+
                          (Protocol) + C1 + C2 + C3 + C4+ C5, 
                        data = covs)
   }
   if (grepl("caudate",rse.name) | grepl("dACC",rse.name) | grepl("sACC",rse.name)) {
      #Skipping Dataset/Protocol correction 
    mod = model.matrix(~ (Dx) + Age +  (Sex) + (mitoRate) + (rRNA_rate) + (totalAssignedGene) + (RIN) +
                          C1 + C2 + C3 + C4+ C5,
                       data = covs)
  }
  if(grepl("amygdala", rse.name)){
     #Skipping Protocol correction for Amygdala
     mod = model.matrix(~  (Dx) + Age +  (Sex) + mitoRate + rRNA_rate + totalAssignedGene + RIN +
                           (Dataset)+ C1 + C2 + C3 + C4+ C5,
                        data = covs)
  }
  
  # ## Fabiana remove constant columns
  if(length(which(colSums(mod)==0))!=0){  mod = mod[,-which(colSums(mod)==0)]
  }
  
 
  ## Compute PCs
  pcobj = prcomp(t(exp),center = T, scale. = T)
  
  ## Check correlation between PCs and neu
  library(corrplot)
  pearson_corr.pc1 = cor.test(as.data.frame(colData(rse))$neu,pcobj$x[,1], method = "pearson")
  pearson_corr.pc2 = cor.test(as.data.frame(colData(rse))$neu,pcobj$x[,2], method = "pearson")
  pearson_corr.pc3 = cor.test(as.data.frame(colData(rse))$neu,pcobj$x[,3], method = "pearson")
 
  corrplot(cor(as.data.frame(colData(rse))$neu, pcobj$x[,1:10]))
  print(paste0('cor: ',pearson_corr.pc1$estimate,' pval: ',pearson_corr.pc1$p.value))
  print(paste0('cor: ',pearson_corr.pc2$estimate,' pval: ',pearson_corr.pc2$p.value))
  print(paste0('cor: ',pearson_corr.pc3$estimate,' pval: ',pearson_corr.pc3$p.value))
  
  
  mod   = cbind(mod, pcobj$x[,seq(from = 1, length.out = 3)])
  
  #Check of multi-colinearity between known covariates and PCs
  print(sum(alias(lm(t(exp)~mod))$Complete))
  
  ## add PCs to colData
  
  stopifnot(identical(rownames(colData(rse)),rownames(pcobj$x[,seq(from = 1, length.out = 3)])))
  colData(rse) = cbind(colData(rse),pcobj$x[,seq(from = 1, length.out = 3)])
  
  cleaned_matrix = cleaningY(exp, mod = mod, P = protect)
  assays(rse)$svacleaned = cleaned_matrix
  
  check.shapiro = apply(cleaned_matrix,1, function(r) shapiro.test(r)$p.value)
  print(paste0("Non-normal genes (p <0.05) before ranknorm: ", sum(check.shapiro <0.05), " / ",length(check.shapiro)))
  
  exp_ranknorm           = t(apply(cleaned_matrix, 1, RNOmni::RankNorm))  #rankNorm each gene in the matrix
  assays(rse)$ranknorm   = exp_ranknorm
  
  exp_zscores             = t(apply(exp_ranknorm, 1, scale)) 
  colnames(exp_zscores)   = colnames(exp_ranknorm)
  assays(rse)$zscores     = exp_zscores
  
  check.shapiro = apply(exp_zscores,1, function(r) shapiro.test(r)$p.value)
  print(paste0("Non-normal genes (p <0.05) after ranknorm: ", sum(check.shapiro <0.05), " / ",length(check.shapiro)))
  
  
  
  return(rse)
}


##Regress confounders and apply genewise rank-normalisation
rse.cleaned.big_logTPM = list()
for (l in names(rse_preprocessed)){
   rse.cleaned.big_logTPM[[l]]    = get.cleaned.obj(rse = rse_preprocessed[[l]], protect=1,
                                                                 rse.name =  l,   assay.to.clean = "logTPM")
}


saveRDS(rse.cleaned.big_logTPM ,"cleaned_data.rds")




