library(purrr)
library(dplyr)
# Set working directory - adjust path as needed
# setwd("path/to/STEP7.Training_INGENE/")


dataset = "your_training_dataset_name"
regions = c("dlpfc")


walk(regions, function(region){
  files = list.files(sprintf("%s/output/summary/%s", getwd(), dataset,region))
  extra = bind_rows(map(files, function(file){
    read.delim(sprintf(sprintf("%s/output/summary/%s/%s", getwd(), dataset,region,file)))
  }))
  save(extra, file=sprintf("%s/output/%s", getwd(), dataset,paste0("INGENE_extraDB_",region,".RData")))
  
  
})
