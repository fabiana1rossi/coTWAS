library(dplyr)
library(readxl)
library(readr)

# Set base directories
# base_dir = "path/to/conditional/results/"
# metadata_path = "path/to/metadata.RData"
base_dir = "./result/conditional.results/"
metadata_path = "./metadata.RData"

check = list.files(file.path(base_dir, "cohort.list"))
check = gsub("\\.rds","",limma::strsplit2(check,"_")[,2])

# Function to calculate root-effective sample size
calc_neff_root <- function(n_case, n_ctrl) {
  neff <- 4 / (1/n_case + 1/n_ctrl)
  sqrt(neff)
}

##### load metadata
load(metadata_path)

##### load Loci
load(file.path(base_dir, "Loci.RData"))
ll = unique(Loci$Region)

purrr::walk(ll[ll %in% check],function(locus){
  print(locus)
  
  #### load data
  cohort.list <- readRDS(paste0(base_dir, "cohort.list/cohort.list_",locus,".rds"))
  stopifnot(identical(names(cohort.list),names(metadata.PRS.sub)))
  stopifnot(sapply(cohort.list,nrow) == sapply(metadata.PRS.sub,nrow))
  
  ##### set weights
  weights = sapply(metadata.PRS.sub,function(x){
    calc_neff_root(n_case = sum(x$Dx == 2),
                   n_ctrl = sum(x$Dx == 1))
  })
  
  ##### create a summary with variance, mean and weights for each gene for each cohort
  cohort.list = Filter(function(x)length(x) > 2,cohort.list)
  summary = sapply(names(cohort.list),function(x){
    df = cohort.list[[x]]
    vv = sapply(as.data.frame(df[,grep("ENSG",colnames(df))]),var)
    mean = sapply(as.data.frame(df[,grep("ENSG",colnames(df))]),mean)
    weight = weights[[x]]
    data.frame(ID = grep("ENSG",colnames(df),value = T),
               var = vv,
               mean = mean,
               weight = weight,
               cohort = x,
               N = nrow(df),row.names = NULL)
  },simplify = F)
  summary = Reduce(rbind,summary)
  summary = summary[complete.cases(summary$var),]
  
  ##### create pooled weighted variance (GVAR)
  summary <- summary %>%
    group_by(ID) %>%
    summarise(
      GVAR = sum(weight * var) / sum(weight),
      N = sum(N)
    )
  
  ##### add Beta, pvalues, N to summary data
  stopifnot(all(summary$ID %in% Loci$tissue.gene))
  summary[,c("beta","SE","p","p.adjusted")] = Loci[match(summary$ID,Loci$tissue.gene),
                                                      c("beta.logistic","se.logistic","meta.logistic.p","meta.p.logistic.adjusted.gwide")]
  summary$af = summary$pos = 0
  
  save(summary, file = paste0(base_dir, "summary/summary.input_",locus,".RData"))
  
})
