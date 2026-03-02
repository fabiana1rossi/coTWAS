# Set base directories
# base_dir = "path/to/conditional/results/"
# metadata_path = "path/to/metadata.RData"
base_dir = "./result/conditional.results/"
metadata_path = "./metadata.RData"

##### load metadata
load(metadata_path)

##### load Loci
load(file.path(base_dir, "Loci.RData"))
to.keep = names(table(Loci$Region))[table(Loci$Region) > 1]

purrr::walk(to.keep,function(locus){
  print(locus)
  #### load data
  cohort.list <- readRDS(paste0(base_dir, "cohort.list/cohort.list_",locus,".rds"))
  stopifnot(identical(names(cohort.list),names(metadata.PRS.sub)))
  stopifnot(sapply(cohort.list,nrow) == sapply(metadata.PRS.sub,nrow))
  
  # Function to calculate root-effective sample size
  calc_neff_root <- function(n_case, n_ctrl) {
    neff <- 4 / (1/n_case + 1/n_ctrl)
    sqrt(neff)
  }
  
  ##### set weights
  weights = sapply(metadata.PRS.sub,function(x){
    calc_neff_root(n_case = sum(x$Dx == 2),
                   n_ctrl = sum(x$Dx == 1))
  })
  
  # Placeholder for weighted sum
  # sum_weighted_corr <- NULL
  
  # Loop over cohorts for genes available in all cohorts
  cohort.list = Filter(function(x)length(x) == 2+(sum(Loci$Region %in% locus)),cohort.list)
  to.get = Reduce(intersect,sapply(cohort.list,function(x)grep("ENSG",colnames(x),value = T),simplify = F))
  
  weighted_correlation = purrr::imap(cohort.list,function(x,y){
    print(y)
    grex <- x
    grex = as.matrix(grex[,to.get])
    grex_cor <- cor(grex, use = "pairwise.complete.obs")
    
    weight <- weights[[y]]
    
    sum_weighted_corr <- grex_cor * weight
  })
  
  # Final weighted average correlation matrix
  sum_weighted_corr = Reduce('+',weighted_correlation)
  avg_correlation_matrix <- sum_weighted_corr / sum(weights)
  
  save(avg_correlation_matrix, 
       file = paste0(base_dir, "corr.matrix/cor.matrix_",locus,".RData"))
  
})




