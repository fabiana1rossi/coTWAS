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
rm(testing_data_lst)

################################################################################
############# COMMON GENES #####################################################
################################################################################
regions = names(GTEX.predictions.ANOVA)
combos = names(GTEX.predictions.ANOVA$dlpfc$common_predictions) #[c(1,2,3,4)]


Common.Predictions.logTPM = sapply(regions, function(region){
  if(region=="amygdala"){combos=combos[!grepl("E",combos)]}
  sapply(combos, function(combo){
    print(paste0(region, ": ", combo))
    All.df = data.frame(GTEX.predictions.ANOVA[[region]]$common_predictions[[combo]])
    gtex.df = data.frame(colData(gtex[[region]]))
    gtex.df = gtex.df[gtex.df$gtex.subjid %in% rownames(All.df),]
    All.df = All.df[rownames(All.df) %in% gtex.df$gtex.subjid,]
    cov.df = gtex.df[match(rownames(All.df), gtex.df$gtex.subjid),c("gtex.subjid", "mean.age", paste0("C", 1:5),
                                                                    "gtex.smrin", "gtex.smrrnart", "gtex.smmaprt","gtex.smtsisch",
                                                                    "neu","gtex.sex" )]
    # cov.df = cov.df[!is.na(cov.df$gtex.subjid),]
    covs = colnames(cov.df)[-1]
    genes = unique(strsplit2(colnames(All.df), split = '\\.')[,2])
    combo.df = sapply(genes, function(gene){
      print(paste0("Gene: ", which(genes %in% gene), "/", length(genes)))
      if(combo %in% c("EMI","EM","EI",
                      "MIC","MC","IC" )){
            
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
        pred.gene = All.df[,c(grep(gene,colnames(All.df)))]
        df.lm = cov.df
        stopifnot(identical(df.lm$gtex.subjid,rownames(pred.gene)))
        
        df.lm$true.expr = true.expr
        df.lm = data.frame(cbind(cov.df,All.df[,c(grep(gene,colnames(All.df)))]),  true.expr)
        MODULE.gene = grep("MODULE",colnames(df.lm),value=T)
        INGENE.gene = grep("INGENE",colnames(df.lm),value=T)
        EpiXcan.gene = grep("EpiXcan",colnames(df.lm),value=T)
        CIS.gene = grep("CIS",colnames(df.lm),value=T)
        
        cis.pred = c(EpiXcan.gene, CIS.gene)
        trans.pred = c(MODULE.gene, INGENE.gene)
        model.cov = paste0(covs, collapse = "+")
        model.cis = paste0(c(covs, cis.pred), collapse = "+")
        model.trans = paste0(c(covs, cis.pred,trans.pred), collapse = "+")
        
        cov.conf = as.formula(paste0("true.expr ~ ",model.cov))
        cov.fit = lm(cov.conf, data = df.lm)
        summary(cov.fit)
        
        cis.conf = as.formula(paste0("true.expr ~ ",model.cis))
        cis.fit = lm(cis.conf, data = df.lm)
        summary(cis.fit)
        
        trans.conf = as.formula(paste0("true.expr ~ ",model.trans))
        trans.fit = lm(trans.conf, data = df.lm)
        summary(trans.fit)
        
        anova.cov = aov(cov.conf, data = df.lm)
        anova.cis = aov(cis.conf, data = df.lm)
        anova.trans = aov(trans.conf,data = df.lm)
        
        stats = list(trans.lm = tidy(trans.fit),
                     cis.lm = tidy(cis.fit),
                     cov.lm = tidy(cov.fit),
                     # combined.fit = trans.fit,
                     trans.aov = tidy(anova.trans),
                     cis.aov = tidy(anova.cis),
                     cov.aov = tidy(anova.cov),
                     trans.lm.stats = glance(trans.fit),
                     cis.lm.stats = glance(cis.fit),
                     cov.lm.stats = glance(cov.fit),
                     trans.aov.stats = glance(anova.trans),
                     cis.aov.stats = glance(anova.cis),
                     cov.aov.stats = glance(anova.cov),
                     trans.on.cis.anova.F = anova(trans.fit, cis.fit)$F[2],
                     trans.on.cis.anova.p = anova(trans.fit, cis.fit)$`Pr(>F)`[2],
                     trans.on.cov.anova.F = anova(trans.fit, cov.fit)$F[2],
                     trans.on.cov.anova.p = anova(trans.fit, cov.fit)$`Pr(>F)`[2],
                     cis.on.cov.anova.F = anova(cis.fit, cov.fit)$F[2],
                     cis.on.cov.anova.p = anova(cis.fit, cov.fit)$`Pr(>F)`[2],
                     adj.r.trans = summary(trans.fit)$adj.r.squared,
                     adj.r.cis = summary(cis.fit)$adj.r.squared,
                     adj.r.cov = summary(cov.fit)$adj.r.squared,
                     Vuong.trans.on.cis = vuongtest(trans.fit, cis.fit)[[2]],
                     Vuong.trans.on.cov = vuongtest(trans.fit, cov.fit)[[2]],
                     Vuong.cis.on.cov = vuongtest(cis.fit, cov.fit)[[2]],
                     AIC.trans.on.cis = icci(trans.fit, cis.fit)[[5]][[1]],
                     AIC.cis.on.trans = icci(trans.fit, cis.fit)[[5]][[2]],
                     AIC.trans.on.cov = icci(trans.fit, cov.fit)[[5]][[1]],
                     AIC.cov.on.trans = icci(trans.fit, cov.fit)[[5]][[2]],
                     AIC.cis.on.cov = icci(cis.fit, cov.fit)[[5]][[1]],
                     AIC.cov.on.cis = icci(cis.fit, cov.fit)[[5]][[2]],
                     Cox.trans.on.cis = coxtest(trans.fit, cis.fit)[[4]][[2]],
                     Cox.trans.on.cov = coxtest(trans.fit, cov.fit)[[4]][[2]],
                     Cox.cis.on.cov = coxtest(cis.fit, cov.fit)[[4]][[2]],
                     jtest.trans.on.cis = jtest(trans.fit, cis.fit)[[4]][[2]],
                     jtest.trans.on.cov = jtest(trans.fit, cov.fit)[[4]][[2]],
                     jtest.cis.on.cov = jtest(cis.fit, cov.fit)[[4]][[2]]
        )
      }else{ # "MI_common_prediction" 
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
        pred.gene = All.df[,c(grep(gene,colnames(All.df)))]
        df.lm = cov.df
        stopifnot(identical(df.lm$gtex.subjid,rownames(pred.gene)))
        
        df.lm$true.expr = true.expr
        df.lm = data.frame(cbind(cov.df,All.df[,c(grep(gene,colnames(All.df)))]),  true.expr)
        MODULE.gene = grep("MODULE",colnames(df.lm),value=T)
        INGENE.gene = grep("INGENE",colnames(df.lm),value=T)
        
        trans.pred = c(MODULE.gene, INGENE.gene)
        model.cov = paste0(covs, collapse = "+")
        model.trans = paste0(c(covs,trans.pred), collapse = "+")
        
        cov.conf = as.formula(paste0("true.expr ~ ",model.cov))
        cov.fit = lm(cov.conf, data = df.lm)
        summary(cov.fit)
        
        trans.conf = as.formula(paste0("true.expr ~ ",model.trans))
        trans.fit = lm(trans.conf, data = df.lm)
        summary(trans.fit)
        
        anova.trans = aov(trans.conf,data = df.lm)
        anova.cov = aov(cov.conf,data = df.lm)
        
        stats = list(trans.lm = tidy(trans.fit),
                     cov.lm = tidy(cov.fit),
                     # combined.fit = trans.fit,
                     trans.aov = tidy(anova.trans),
                     cov.lm = tidy(cov.fit),
                     cov.aov = tidy(anova.cov),
                     trans.lm.stats = glance(trans.fit),
                     trans.aov.stats = glance(anova.trans),
                     cov.lm.stats = glance(cov.fit),
                     cov.aov.stats = glance(anova.cov),
                     trans.on.cov.anova.F = anova(trans.fit, cov.fit)$F[2],
                     trans.on.cov.anova.p = anova(trans.fit, cov.fit)$`Pr(>F)`[2],
                     adj.r.trans = summary(trans.fit)$adj.r.squared,
                     adj.r.cov = summary(cov.fit)$adj.r.squared,
                     Vuong.trans.on.cov = vuongtest(trans.fit, cov.fit)[[2]],
                     AIC.trans.on.cov = icci(trans.fit, cov.fit)[[5]][[1]],
                     AIC.cov.on.trans = icci(trans.fit, cov.fit)[[5]][[2]],
                     Cox.trans.on.cov = coxtest(trans.fit, cov.fit)[[4]][[2]],
                     jtest.trans.on.cov = jtest(trans.fit, cov.fit)[[4]][[2]]
        )
      }
    },simplify = F, USE.NAMES=T) 
    save(combo.df, file = paste0("predictions/LIBD_training/GTeX.v9/Combined/MLE_common_",region,"_",combo, ".RData"))
    return(combo.df)
  },simplify = F, USE.NAMES=T)
},simplify = F, USE.NAMES=T)

save(Common.Predictions.logTPM, file = "predictions/LIBD_training/GTeX.v9/Combined/MLE_common.RData")
