library(dplyr)
library(purrr)

# Set base directory for conditional analysis results
# base_dir = "path/to/conditional/results/"
base_dir = "./result/conditional.results/"

##### load Loci
load(file.path(base_dir, "Loci.RData"))
filter1 = names(table(Loci$Region))[table(Loci$Region) == 1]
filter2 = gsub("\\.RData","",limma::strsplit2(grep("results_",list.files(base_dir),value = T),"_")[,2])

#########
Loci_1 = purrr::map_dfr(filter1,function(x){
  load(paste0(base_dir, "summary/summary.input_",x,".RData"))
  colnames(summary)[1] = "rsid"
  summary$pos = summary$af = NULL
  summary[, c("beta_old", "se_old", "p_old", "beta_new",  "se_new", "p_new")] = cbind(summary[,c("beta","SE","p")],summary[,c("beta","SE","p")])
  summary$z_old = summary$z_new = summary$beta/summary$SE
  summary$Adjusted_Beta = summary$beta_new * sqrt(summary$GVAR)
  summary$Adjusted_OR <- exp(summary$Adjusted_Beta)
  return(summary)
})

Loci_2 = purrr::map_dfr(filter2,function(x){
  load(paste0(base_dir, "results_",x,".RData"))
  gene_stats_cond$pos.x = gene_stats_cond$pos.y = gene_stats_cond$af = NULL
  return(gene_stats_cond)
})

Loci_cond = rbind(Loci_2,Loci_1[,colnames(Loci_2)])
Loci_cond$p_new.adjusted = p.adjust(Loci_cond$p_new, method = "BH")
Loci_cond$Region = Loci$Region[match(Loci_cond$rsid,Loci$tissue.gene)]
save(Loci_cond, file = file.path(base_dir, "Loci_cond.RData"))

###### set cohort list
# Load metadata from parent directory
# metadata_path = "path/to/metadata.RData"
metadata_path = "./metadata.RData"
load(metadata_path)
load(file.path(base_dir, "Loci_cond.RData"))
Loci_cond = Loci_cond[Loci_cond$p_new.adjusted <= .01,]
load(file.path(base_dir, "Loci.RData"))
Loci_cond$Region = Loci$Region[match(Loci_cond$rsid,Loci$tissue.gene)]
rm(Loci,metadata)

cohort_list = sapply(unique(Loci_cond$Region), function(locus){
  
  ll = readRDS(paste0(base_dir, "cohort.list/cohort.list_",locus,".rds"))
  ll = lapply(ll,function(x) x[,colnames(x) %in% c("FID","IID",Loci_cond$rsid[Loci_cond$Region %in% locus])])
  
},simplify = F)

merged_data_list <- names(cohort_list[[1]]) %>%
  map(~ {
    # .x here refers to the current dataframe name (e.g., "gene_data")
    print(.x)
    # Create a list of the corresponding dataframes from each main list
    list_of_corresponding_dfs <- map(cohort_list, `[[`, .x)
    
    # Iteratively merge all dataframes in the list
    reduce(
      list_of_corresponding_dfs,
      full_join,
      by = c("FID","IID")
    )
  }) %>%
  set_names(names(cohort_list[[1]]))

stopifnot( all(sapply(merged_data_list,nrow) == sapply(metadata.PRS.sub,nrow)))

merged_data_list = purrr::map2(merged_data_list,metadata.PRS.sub, function(x,y){  merge(x,y,by=c("FID","IID"))  })
save(merged_data_list, file = file.path(base_dir, "cohort_cond.RData"))




