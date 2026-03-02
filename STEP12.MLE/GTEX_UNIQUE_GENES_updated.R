# Set working directory to project root
# drive = "path/to/project/root"
# setwd(drive)
base_dir = getwd()
load("./predictions/LIBD_training/GTeX.v9/Combined/predictions_testing_data_lst.RData")
gtex = readRDS("../../dataset/GTEx/expression/rse_GTEX_cleaned_gtex.rds")
gtex$sACC = gtex$ACC
gtex$dACC = gtex$ACC

library(tidyr)
library(dplyr)
library(limma)
library(broom)
library(nonnest2) 
library(lmtest)
library(caret)
library(SummarizedExperiment)

GTEX.predictions.ANOVA = testing_data_lst
models = names(GTEX.predictions.ANOVA$dlpfc$unique_predictions)
regions = names(GTEX.predictions.ANOVA)

################################################################################
############# UNIQUE GENES #####################################################
################################################################################


Unique.Predictions.logTPM = sapply(regions, function(region){
  sapply(models, function(model){
    print(paste0(region, ": ", model))
    All.df = data.frame(GTEX.predictions.ANOVA[[region]]$unique_predictions[[model]])
   
    gtex.df = data.frame(colData(gtex[[region]]))
    gtex.df = gtex.df[gtex.df$gtex.subjid %in% rownames(All.df),]
    All.df = All.df[rownames(All.df) %in% gtex.df$gtex.subjid,]
    cov.df = gtex.df[match(rownames(All.df), gtex.df$gtex.subjid),c("gtex.subjid", "mean.age", paste0("C", 1:5),
                                                                    "gtex.smrin", "gtex.smrrnart", "gtex.smmaprt","gtex.smtsisch",
                                                                    "neu","gtex.sex" )]
    covs = colnames(cov.df)[-1]
    genes = unique(colnames(All.df))
   
    print(paste0("Number of Genes: ", length(genes)))
    sapply(genes, function(gene){
      print(paste0("Gene No. ", which(genes %in% gene), "/ ", length(genes)))
      
      predict.df = as.data.frame(t(as.data.frame(gtex[[region]]@assays@data$logTPM)))
      predict.df = predict.df[which(rownames(predict.df) %in% gtex.df$external_id),]
      gtex.df    = gtex.df[which(gtex.df$external_id %in% rownames(predict.df)),]
      
      predict.df = predict.df[match(gtex.df$external_id,rownames(predict.df)),]
      stopifnot(identical(gtex.df$external_id,rownames(predict.df)))
      rownames(predict.df) = gtex.df$gtex.subjid
      
      predict.df = predict.df[match(rownames(All.df),rownames(predict.df)),]
      stopifnot(identical(rownames(All.df),rownames(predict.df)))
      colnames(predict.df) = strsplit2(colnames(predict.df),"\\.")[,1]
      
      true.expr = predict.df[[gene]]
      
      pred.gene = All.df[,c(grep(gene,colnames(All.df))),drop=F]
      df.lm = cov.df
      stopifnot(identical(df.lm$gtex.subjid,rownames(pred.gene)))
      
      df.lm = data.frame(cbind(pred.gene[,1], cov.df, true.expr))
      colnames(df.lm)[1] = gene
      model.conf = paste0(c(covs, gene), collapse = "+")
      model.fm = as.formula(paste0("true.expr ~ ",model.conf))
      model.fit = lm(model.fm, data = df.lm)
      summary(model.fit)
      
      cov.conf = paste0(c(covs), collapse = "+")
      cov.fm = as.formula(paste0("true.expr ~ ",cov.conf))
      cov.fit = lm(cov.fm, data = df.lm)
      summary(cov.fit)
      
      
      anova.model = aov(model.fm, data = df.lm)
      anova.cov = aov(cov.fm,data = df.lm)
      
      stats = list(model.lm = tidy(model.fit),
                   cov.lm = tidy(cov.fit),
                   model.lm.stats = glance(model.fit),
                   cov.lm.stats = glance(cov.fit),
                   model.aov = tidy(anova.model),
                   cov.aov = tidy(anova.cov),
                   model.aov.stats = glance(anova.model),
                   cov.aov.stats = glance(anova.cov),
                   anova.F = anova(model.fit, cov.fit)$F[2],
                   anova.p = anova(model.fit, cov.fit)$`Pr(>F)`[2],
                   adj.r.model = summary(model.fit)$adj.r.squared,
                   adj.r.cov = summary(cov.fit)$adj.r.squared,
                   adj.r.diff = summary(model.fit)$adj.r.squared - summary(cov.fit)$adj.r.squared,
                   Vuong = vuongtest(model.fit, cov.fit)[[2]],
                   AIC.trans = icci(model.fit, cov.fit)[[5]][[1]],
                   AIC.cis = icci(model.fit, cov.fit)[[5]][[2]],
                   Cox = coxtest(model.fit, cov.fit)[[4]][[2]],
                   jtest = jtest(model.fit, cov.fit)[[4]][[2]]
      )
      
    },simplify = F, USE.NAMES=T) 
  },simplify = F, USE.NAMES=T)
},simplify = F, USE.NAMES=T)

save(Unique.Predictions.logTPM, file = "predictions/LIBD_training/GTeX.v9/Combined/MLE_unique.RData")
