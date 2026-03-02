library(purrr)
library(dplyr)
library(RSQLite)

# Get base directory for relative paths
# base_dir = "path/to/project/root"
base_dir = getwd()


# Configuration settings
CONFIG    = list(
    train_name    = "your_training_dataset_name",
    db_file_path  = file.path(base_dir,"./PredictDB_LD/MODULE/"),
    regions       = c("dlpfc"),
    out_dir       = file.path(base_dir, "PredictDB_LD/MODULE/")
  )





extractDBinfo = function(brain.region, train_name, db_file_path, out_dir){
   print(brain.region)
  
  ## get DB
  region_db_file_path = file.path(db_file_path, brain.region)
  dbs                 = list.files(region_db_file_path)
 # dbs                 = dbs[grepl(paste0(brain.region,"_",train_name),dbs)]
  
  net_lst = map(dbs, ~{
    print(.x)
    name_drv_file = file.path(region_db_file_path, .x)
    
    ## *db file is a SQLite database.
    driver = dbDriver('SQLite')
    conn   = dbConnect(drv = driver,  name_drv_file)
    
    ## Extract weight Table  
    extra   = dbReadTable(conn, "extra")
    weights = dbReadTable(conn, "weights")
    
    #weights = weights %>% mutate(network = unique(extra$network))
    dbDisconnect(conn)
    return(list(extra=extra,weights=weights))
    
  })
  
  
  ## Concatenate df for brain region
  extras   = do.call(rbind,map(net_lst, ~ .x[["extra"]]))
  weights  = do.call(rbind,map(net_lst, ~ .x[["weights"]]))
  
  ## Save results
  save(extras,weights, file = paste0(out_dir, "/MODULE.",brain.region,".extras_weights.RData"))  
  
}



for(region in CONFIG$regions){
  extractDBinfo(region, CONFIG$train_name, CONFIG$db_file_path, CONFIG$out_dir)
}
