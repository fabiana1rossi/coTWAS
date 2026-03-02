###############################################################################################
# Script to Preprocess tissue rse data and remove outliers
# This script performs:
# 1. Loads RSE (SummarizedExperiment) data and processes expression assays
# 2. Filters samples by age, race, diagnosis, and genomic information
# 3. Removes outlier samples using IAC (Inter-Array Correlation) analysis
# 4. Removes outlier genes based on expression thresholds
# 5. Converts RPKM to TPM and calculates cell type proportions
# This script is adapted from: https://github.com/LieberInstitute/Brain_WGCNA/tree/main/preprocess
###############################################################################################

options(java.parameters = "-Xmx8000m")  # Required for xlsx package
library(SummarizedExperiment)           # check this link: https://f1000research.com/articles/6-1558/v1 , https://www.bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
library(recount)                        # check this paper https://www.nature.com/articles/nbt.3838
# library(jaffelab)                       # vignette http://127.0.0.1:13522/library/recount/doc/recount-quickstart.html
library(limma)
library(dendextend)
library(purrr)
library(dplyr)
#library(preprocessCore)
library(pheatmap)
library(RNOmni)
#library(quantro)
library(BRETIGEA)

set.seed(123)

setwd('')


#Read genomic eigenvariables (EUR)
load("./GE.RData") ## your GE data
gnm_eur = GE
colnames(gnm_eur)[2] = 'IID'
gnm_eur = gnm_eur[!duplicated(gnm_eur$IID),]
gnm = gnm_eur
gnm = gnm[,1:7] 



getRSEdata = function(path){
  rse = list()
  load(path)
  rse$dlpfc = subset(rse_gene,colData(rse_gene)$Region == 'DLPFC',subset = TRUE)
  rse$caudate <- subset(rse_gene, colData(rse_gene)$Region == "Caudate", subset=TRUE)
  rse$hippo <- subset(rse_gene, colData(rse_gene)$Region == "HIPPO",subset = TRUE)
  rse$amygdala <- subset(rse_gene, colData(rse_gene)$Region == "Amygdala",subset = TRUE)   
  rse$dACC       = subset(rse_gene, colData(rse_gene)$Region == "dACC",subset = TRUE)
  rse$sACC       = subset(rse_gene, colData(rse_gene)$Region == "sACC",subset = TRUE)      
 
  return(rse)
}

#Load all RSE data in a list
rse_raw_list = getRSEdata(path="./merged_expression_assays.rda") ##  your rda data to process

#Save gene_map/gene_info for later reference:
gene_map           = as.data.frame(rowData(rse_raw_list$dlpfc), stringsAsFactors = F)
gene_map           = gene_map[!duplicated(gene_map$ensemblID),]
colnames(gene_map)[grep("ensembl",colnames(gene_map))] = "ensembl"
rownames(gene_map) = gene_map$ensembl


#Combine Genomic Eigen variable information into the rse -colData:
rse_raw_list_updated = map(rse_raw_list,~{
  yy = gnm[match(.$BrNum, gnm$IID), c("PC1","PC2","PC3","PC4","PC5")]
  colnames(yy) = c("C1","C2","C3","C4","C5")
  colData(.) = cbind(colData(.),yy)
  return(.)
})


#Pre-process rse-data
rse_merged_list = rse_raw_list_updated %>% purrr::imap(~{
  #rse_merged = jaffelab::merge_rse_metrics(.) 
  assays(.)$RPKM     = recount::getRPKM(., length_var = "Length", mapped_var = "numMapped")
  assays(.)$logRPKM  = log2(assays(.)$RPKM+1)
  return(.)
})


#Subset RSE objects for desired Samples:
#minAge = min(rse_merged_list$dlpfc@colData$Age);
minAge = 17; maxAge =  max(rse_merged_list$dlpfc@colData$Age);
minRin = 6; diagnosis = c("Control", "SCZD","MDD","Bipolar");  
race = c("CAUC"); 

rse_preprocess_step1 = rse_merged_list %>% purrr::imap(~{
  rse_merged_small = SummarizedExperiment::subset(
    x      = .                                     ,
    subset = TRUE                                  ,                         #----"subset" argument subsets from rowData(rse) i.e. genes  ##(default = T, all rows)
    select =    
                                                                             #----"select" argument subsets from colData(rse) i.e. samples ##(default = T, all columns)
      (Age                >=   minAge       )   &
      (Age                <=   maxAge       )   &
      (Race              %in%  race         )   &
      (Dx                %in%  diagnosis    )   &
      (!is.na(C1)                          )                                #Keep if C1 genomic information is not-missing
  )
  # hist(.@colData@listData[["Age"]], xlab = "Age",main='', col= 'purple')
  return(rse_merged_small)
})



#Outlier sample removal (from each subgroup)
remove_outlier_samples = function(rse, name.rse , sdout = 3){
  exp                 = assays(rse,withDimnames = F)$logRPKM
  IAC                 = as.matrix(dist(t(exp)))            #Sample-sample pairwise euclidean distance/disimilarity (based on all genes)
  
  scaled.IAC          = scale(matrixStats::rowMeans2(IAC))
  names(scaled.IAC)   = rownames(IAC)
  samp.outlier        = abs(scaled.IAC) > sdout
  names(samp.outlier) = rownames(IAC)
  
  IACnew              = IAC[!samp.outlier,!samp.outlier]   #After removing outliers
  
  dend = IAC %>% as.dist %>% hclust (method = "average") %>% as.dendrogram
  
  cll = rep("black", length(labels(dend))); names(cll) = labels(dend)
  cll[labels(dend) %in% names(samp.outlier)[samp.outlier]] = "red"
  dend = dend %>% set("labels_col", cll)
  dend = dend %>% set("labels_cex", 0.6)
  
  dendnew = IACnew %>% as.dist %>% hclust (method = "average") %>% as.dendrogram
  dendnew = dendnew %>% set("labels_cex", 0.6)
  
  par(mfrow = c(2,2), mar = c(2.5,4,1.5,1),oma = c(1,1,0.5,0.5))
  y.lim = c(0,max(attributes(dend)$height,attributes(dendnew)$height))
  plot(dend    ,ylim = y.lim)
  plot(dendnew ,ylim = y.lim)
  plot(scaled.IAC[match(labels(dend),names(scaled.IAC))], ylim = c(-6,6), pch = 19, col = cll, ylab = "z-scored sample")
  abline(h=c(sdout,-sdout), col = "red")
  title(name.rse,outer = T,line = -1)
  
  return(samp.outlier)
}

# Remove outlier samples from age-parsed data using logTPM data
rse_outlier_removed = c(unlist(rse_preprocess_step1)) %>% purrr::imap( ~{
  found.outliers = remove_outlier_samples(.x,.y,sdout = 3)
  rse_merged = SummarizedExperiment::subset(
    x      = .x                                     ,
    select = !found.outliers                        ,
    subset = TRUE
  )
  return(list(rse = rse_merged, removed = names(found.outliers)[found.outliers]))    
})


rse_outlier_removed_age_grps = map(rse_outlier_removed,"rse")
samples_removed              = map(rse_outlier_removed, "removed")


#Re-combining age-parsed groups after groupwise outlier removal
rse_outlier_removed_big = list(
  dlpfc   = do.call(cbind,rse_outlier_removed_age_grps[grep("dlpfc"  , names(rse_outlier_removed_age_grps))]),
  hippo   = do.call(cbind,rse_outlier_removed_age_grps[grep("hippo"  , names(rse_outlier_removed_age_grps))]),
  caudate = do.call(cbind,rse_outlier_removed_age_grps[grep("caudate", names(rse_outlier_removed_age_grps))]),
  amygdala = do.call(cbind,rse_outlier_removed_age_grps[grep("amygdala", names(rse_outlier_removed_age_grps))]),
  sACC = do.call(cbind,rse_outlier_removed_age_grps[grep("sACC", names(rse_outlier_removed_age_grps))]),
  dACC = do.call(cbind,rse_outlier_removed_age_grps[grep("dACC", names(rse_outlier_removed_age_grps))])
)


rse_outlier_removed_age_grps = NULL
############


#Outlier gene removal
rse_preprocess_step2 = imap(rse_outlier_removed_big, ~{
  expMat            = assays(.)$RPKM
  median.gene.exprs = matrixStats::rowMedians(expMat)
  names(median.gene.exprs) = rownames(expMat)
  
  #Gene filters
  gene_filter0 = as.vector(rowSums(expMat ==0) <= ncol(expMat)*0.2)            #Genes to keep (not more than 20% samples are zero)
  gene_filter1 = as.vector(median.gene.exprs   >= .1)                          #Genes to keep (median exp > 0.1)
  
  median.gene.exprs_filtered = median.gene.exprs[gene_filter0 & gene_filter1]
  gene_filter2 = as.vector(abs(scale(median.gene.exprs_filtered)) <  3)        #Genes to keep (z-scored_median_exp < 3)
  
  
  genes_to_keep  = names(median.gene.exprs_filtered)[gene_filter2]
  genes_to_keep  = limma::strsplit2(genes_to_keep,'\\|')[,2]
  
  print(rbind(table(gene_filter0), table(gene_filter1), table(gene_filter2)))
  print(paste0(length(genes_to_keep),"/",dim(expMat)[1]))
  
  rse_merged2 =SummarizedExperiment::subset(
     x      = .                                     ,
     subset = (gencodeID  %in% genes_to_keep)     &             
        (!grepl("^MT-",Symbol)        )                                   #Gene symbol starting with 'MT-' (for mitochondrial genes)
  )
   return(rse_merged2)
 })



### convert in TPM
rse_preprocess_step2  = map(rse_preprocess_step2, ~{
  rpkm = .x@assays@data@listData[["RPKM"]]
  assays(.x)$TPM = apply(rpkm, 2, function(x) {(x/sum(x)) * 10^6})
  assays(.x)$logTPM  = log2(assays(.x)$TPM+1)
  .x
})



# Add cell proportion estimate with package BRETIGEA
get.cells.proportions = function(rse,rse.name, assay.name, method, nmarkers){
  exp = assays(rse)[[assay.name]]       #Genes in rows
  rownames(exp) = rowData(rse)$Symbol   #Convert ENSEMBL names to GeneSymbol
  zz1 = BRETIGEA::brainCells(exp, nMarker=nmarkers, species="human", method = method, scale = T) #Scaled estimates
  zz2 = BRETIGEA::brainCells(exp, nMarker=nmarkers, species="human", method = method, scale = F) #Unscaled estimates
  colnames(zz2) = paste0(colnames(zz2),"_unscaled")
  zz = cbind(zz1,zz2)
  return(zz)
}

rse_preprocess_step3 = rse_preprocess_step2 %>% purrr::imap(~{
  nmarkers = 50
  cell.prop.info_SVD = get.cells.proportions(.x, .y, assay.name = "TPM", method = "SVD", nmarkers)
  colData(.x) = cbind(colData(.x), cell.prop.info_SVD)
  return(.x)
})



saveRDS(rse_preprocess_step3, "outlierRemoved.rds")         #List of rse object preprocessed (sample and gene outlier removed. Quantile)



 
