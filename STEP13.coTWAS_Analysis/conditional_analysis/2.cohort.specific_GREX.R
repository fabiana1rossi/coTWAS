require(data.table)
require(readxl)

# Set base directories
# base_dir = "path/to/conditional/results/"
# predictions_dir = "path/to/predictions/"
base_dir = "./result/conditional.results/"
predictions_dir = "./result/predictions/"

##### load Loci
load(file.path(base_dir, "Loci.RData"))

##### set relative cohorts GREX
cohorts = list.files(predictions_dir)
cohorts = unique(limma::strsplit2(cohorts,"_")[,1])
locus = "MHC"
results = Loci[Loci$Region %in% locus,]
cohort.list = purrr::map(unique(results$tissue),function(x){
  
  df = results[results$tissue %in% x,]
  c.list = sapply(cohorts,function(cc){
    
    dd = readRDS(paste0(predictions_dir, cc,"_",x,"_predictions.rds"))
    dd = dd[,colnames(dd) %in% c("FID","IID",df$gene)]
    
  },simplify = F)
})
names(cohort.list) = unique(results$tissue)
cohort.list = purrr::imap(cohort.list,function(x,y){
  
  lapply(x,function(z){ colnames(z)[grep("ENSG",colnames(z))] = paste0(colnames(z)[grep("ENSG",colnames(z))],".",y); z })
  
})
cohort.list = purrr::transpose(cohort.list)
cohort.list = lapply(cohort.list,function(x) Reduce(cbind,x))
cohort.list = lapply(cohort.list,function(x){ x[,!duplicated(colnames(x))] })
saveRDS(cohort.list,
        file = paste0(base_dir, "cohort.list/cohort.list_",locus,".rds"))



# ##### explore GREX across cohorts
# cohort.list <- readRDS(file.path(base_dir, "cohort.list.rds"))
# common.genes = Reduce(intersect,lapply(cohort.list,colnames))
# common.genes = common.genes[!common.genes %in% c("FID","IID")]
# cohort.list.sub = lapply(cohort.list,function(x) x[,match(common.genes,colnames(x))])
# cohort.list.sub = purrr::transpose(cohort.list.sub)
# pdf(file.path(base_dir, "boxplot.pdf"),width = 14,height = 10)
# par(mar=c(1,1,1,1),mfrow = c(6,1))
# purrr::walk(cohort.list.sub,function(x) boxplot(x))
# dev.off()
# 
# cohort.var = sapply(cohort.list.sub,function(x) sapply(x,var))
# cohort.median = sapply(cohort.list.sub,function(x) sapply(x,median))
# cohort.var.plot = reshape2::melt(cohort.var)
# cohort.median.plot = reshape2::melt(cohort.median)
# 
# require(ggpubr)
# p.var = ggdensity(cohort.var.plot[cohort.var.plot$Var2 %in% common.genes[1:100],], x = "value", fill = "lightgray", rug = F) +
#   facet_wrap(.~Var2, scales = "free")
# p.var
# 
# p.median = ggdensity(cohort.median.plot[cohort.median.plot$Var2 %in% common.genes[1:100],], x = "value", fill = "lightgray", rug = F) +
#   facet_wrap(.~Var2, scales = "free")
# p.median
# 
# pdf(file.path(base_dir, "density.first100.pdf"),width = 16,height = 12)
# print(p.var)
# print(p.median)
# dev.off()
# 
# ################# corrplot in each cohort
# library(heatmap3)
# cohort.list <- readRDS(file.path(base_dir, "cohort.list.rds"))
# common.genes = Reduce(intersect,lapply(cohort.list,colnames))
# common.genes = common.genes[!common.genes %in% c("FID","IID")]
# cohort.list.sub = lapply(cohort.list,function(x) x[,match(common.genes,colnames(x))])
# 
# pdf(file.path(base_dir, "corrplot.pdf"),width = 16,height = 12)
# purrr::iwalk(cohort.list.sub,function(x,y) { 
#   heatmap3(cor(x),useRaster=TRUE, main = y)
# })
# dev.off()





