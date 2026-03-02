
require(rcompanion)
require(broom)
require(dplyr)
require(mgsub)
# library(pbapply)


# Set working directory to project root
# setwd("path/to/project/root")


tissues = c("dlpfc","dACC","sACC")
sites   = c("CMC_EUR")
out.dir = c("./predictions/Combined/")


## load genes
load("./predictions/Combined/EMIC_Common_Unique_updated.RData")
genes = purrr::map2(full$common,full$unique,function(x,y) rbind(x,y))

MIEC_Combined_Predictions_lst = map(sites, function(site){
  region.combined_df = map(tissues %>% set_names(.,.), function(tissue){
    
  print(paste0(site,": ", toupper(tissue)))
  
  # ## load predictions
  predictions.common.unique  = get(load(paste0(out.dir,"/predictions_testing_data_lst.RData")))
  predictions                = predictions.common.unique[[tissue]]
  
  load(paste0("./predictions/Combined/Common_Predictions_GTEx_",tissue,"_replicable.RData"))
  model = lapply(model,purrr::transpose)
 
  #### run analysis 
  Combined.Common_genes      = sapply(genes[[tissue]]$gene,function(gene){
    print(gene)
    
    # print(gene)
    #### check matching between GTEx and site
    mm.name = genes[[tissue]]$model[genes[[tissue]]$gene %in% gene]
    models = limma::strsplit2(mm.name,"")[1,]

    predictions_models = predictions$common_predictions[[mm.name]]
    if(is.null(predictions_models)){message("NULL");return(rep(NA, nrow( predictions$common_predictions$EMI$MODULE)))}
    names(predictions_models)[grepl("MODULE",names(predictions_models))] = "M"
    names(predictions_models)[grepl("INGENE",names(predictions_models))] = "I"
    names(predictions_models)[grepl("EpiXcan",names(predictions_models))] = "E"
    names(predictions_models)[grepl("CIS",names(predictions_models))] = "C"
    
    df = as.data.frame(matrix(data = NA, ncol = length(models)))
    
    df = data.frame(sapply(models, function(x) {
      if (length(which(colnames(predictions_models[[x]]) == gene)) > 0) {
        return(predictions_models[[x]][[gene]])
      } else {
        return(rep(NA, nrow(predictions_models[[x]])))
      }
    }, USE.NAMES = FALSE, simplify = TRUE))
    
    colnames(df) = models
    df = df[, apply(df, 2, function(x) !all(is.na(x))),drop=F]
    
    if(ncol(df)==0){return(as.data.frame(matrix(data = NA, nrow = nrow(df))))}
    
    models.available = colnames(df)
    check.m          = paste0(models.available,collapse = "") == mm.name
    check.length.m   = length(models.available) > 1

    
    #### run for different conditions
    if(check.m & check.length.m) {

        ### predict model from GTEx
        colnames(df) = mgsub::mgsub(colnames(df),
                                    c("^E$","^M$","^I$","^C$"),
                                    c(paste0("EpiXcan.",gene),
                                      paste0("MODULE.",gene),
                                      paste0("INGENE.",gene),
                                      paste0("CIS.",gene)))
        mm           = model[[mm.name]][[paste0("model.",mm.name)]][[gene]]
        predictors   = colnames(mm$model)[-1]
        df           = df[,match(predictors,colnames(df))]
        df$predictor = predict(mm,df[,predictors])

      } else if(check.m & !check.length.m){
        predictors = mm.name
        df$predictor = df[[predictors]]

      } else if(!check.m & check.length.m){

        ### predict model from GTEx
        colnames(df) = mgsub::mgsub(colnames(df),
                                    c("^E$","^M$","^I$","^C$"),
                                    c(paste0("EpiXcan.",gene),
                                      paste0("MODULE.",gene),
                                      paste0("INGENE.",gene),
                                      paste0("CIS.",gene)))
        mm.name.available = paste0(models.available,collapse = "")
        mm                = model[[mm.name]][[paste0("model.",mm.name.available)]][[gene]]
        predictors        = colnames(mm$model)[-1]
        df                = df[,match(predictors,colnames(df))]
        df$predictor      = predict(mm,df[,predictors])

      } else if(!check.m & !check.length.m){
        predictors   = paste0(models.available,collapse = "")
        df$predictor = df[[predictors]]
        
      }
    
    ## return gene vector 
    df.gene            = df[,"predictor",drop=F]
    colnames(df.gene) = gene

    return(df.gene)
  },USE.NAMES = FALSE, simplify = TRUE)
  
  Combined.Common_genes           = do.call("cbind",Combined.Common_genes)
  rownames(Combined.Common_genes) = rownames(predictions$common_predictions$EMI$MODULE)
  combined_pred_df                =  as.data.frame(Combined.Common_genes[, apply(Combined.Common_genes, 2, function(x) !all(is.na(x))),drop=F])
  combined_pred_df$FID = combined_pred_df$IID = rownames(combined_pred_df)
  write.table(combined_pred_df, file=paste0(out.dir,tissue,"_MIEC_predictions.txt"), sep = "\t",quote=F)
  
  return(combined_pred_df)
  })
  
  saveRDS(region.combined_df, file=paste0(out.dir,"CMC_EUR_allRegions_MIEC_predictions.rds"))
})



