library(purrr)

load("./metadata.RData")

 tissues = c("dlpfc","hippo","caudate","amygdala","dACC","sACC")
 #tissues = c("sACC")

#### run jobs
purrr::walk(tissues,function(tt){
  arrays = paste0("1-",length(metadata.PRS.sub))
  #arrays = 1  
  print(paste0("running job for: ",tt))
  cmd = paste0("sbatch -a ",arrays,"./2.PGC_TWAS_Analyses.sh ",tt)
  system(cmd)
})




