#' Create SQLite database for prediction models
#' @param brain.region Brain region name
#' @param file_input_summ Directory containing summary files
#' @param file_input_weights Directory containing weights files
#' @param drv_file_path Directory for database files
#' @param net.names Vector of network names
#' @param no_cores Number of cores for parallel processing
#' @param dataset_training Dataset identifier
#' @return List of created database file paths
makeModelDB1 = function(brain.region, file_input_summ, file_input_weights, 
                       drv_file_path, net.names, no_cores, dataset_training) {
  
  make_models = map(net.names %>% set_names(.,.), ~ {
    
    # Read summary file for network
    summary = read.table(file.path(file_input_summ, brain.region, paste0(.x, ".extra.txt")), header = TRUE, stringsAsFactors = FALSE)
    weights = read.table(file.path(file_input_weights, brain.region, paste0(.x, ".weights.txt")), header = TRUE, stringsAsFactors = FALSE)

    # Create directory for database files if it doesn't exist
    saved_db_file_path = file.path(drv_file_path, brain.region)
    if (!dir.exists(saved_db_file_path)) {
      dir.create(saved_db_file_path, recursive = TRUE)
    }
    ## Assign name to DB file
    db_file = file.path(saved_db_file_path, 
                       paste0(dataset_training, 
                             '_', .x, '.db'))
    
    if(file.exists(db_file)) {
      return(db_file)
    }
    
    ## Create DB
    res = makeModelDB2(
      brain.region       = brain.region,
      summary            = summary,
      weights            = weights,
      saved_db_file_path = saved_db_file_path,
      output_drv_name    = basename(db_file),
      no_cores           = no_cores
    )
    
    return(res)
  })
  
  print(make_models)
  return(make_models)
}

#' Process model summaries and weights to create SQLite database
#' @param modules Vector of module names
#' @param brain.region Brain region name
#' @param ethnicity Ethnicity identifier
#' @param file_input_summ Directory containing summary files
#' @param file_input_weights Directory containing weights files
#' @param net.name Network name
#' @param drv_file_path Directory for database files
#' @param output_drv_name Output database filename
#' @param saved_db_file_path Directory for saved database files
#' @param rsID_df Optional rsID mapping dataframe
#' @param no_cores Number of cores for parallel processing
#' @return Path to created database file
makeModelDB2 = function(brain.region, summary, weights, 
                       saved_db_file_path, output_drv_name,  no_cores) {
  #########################################
  # Filter model with significant results
  #########################################
  
  model_summaries = summary %>% dplyr::rename(genename = gene_name, n.snps.in.model = n_snps_in_model, 
                  gene = gene_id)
  model_summaries$pred.perf.qval = model_summaries$pred.perf.R2 = model_summaries$pred.perf.pval = NA 
  
  model_summaries = model_summaries %>% filter(adj_rsq_gene_all_data  >= 0.01 & pval_gene_all_data < .05)

    
   ## Weights Table
  weights = weights %>% dplyr::filter(gene %in% model_summaries$gene) %>%
    dplyr::rename(rsid = snp,ref_allele = ref_vcf, eff_allele = alt_vcf, weight = weights) %>%
    dplyr::select(gene, rsid, varID, ref_allele, eff_allele, weight)
  
  
  ## Check
  weights = weights[!duplicated(weights), ]
  weights$weight = as.numeric(weights$weight)
  weights = na.omit(weights)
 
  ## Filter weights and model summaries to have the same genes
  weights         = weights %>% filter(gene %in% model_summaries$gene)
  model_summaries = model_summaries %>% filter(gene %in% weights$gene)
  
  
  ################################################
  #  CREATE DB ##
  ################################################
  driver = dbDriver("SQLite")
  
  # Connect to a new or existing SQLite database
  con = dbConnect(RSQLite::SQLite(), file.path(saved_db_file_path, output_drv_name))
  
  # Check if tables exist and drop them if they do
  existing_tables = dbListTables(con)
  if("extra" %in% existing_tables) {
    dbExecute(con, "DROP TABLE extra")
  }
  if("weights" %in% existing_tables) {
    dbExecute(con, "DROP TABLE weights")
  }
  
  # Create the extra table
  dbExecute(con, "
  CREATE TABLE extra (
    extra_id INTEGER PRIMARY KEY,
    description TEXT NOT NULL,
    details TEXT
  );
 ")
  
  # Create the weights table
  dbExecute(con, "
  CREATE TABLE weights (
    weight_id INTEGER PRIMARY KEY,
    value REAL NOT NULL,
    unit TEXT NOT NULL,
    extra_id INTEGER,
    FOREIGN KEY(extra_id) REFERENCES extra(extra_id)
  );
 ")
  
 
  ## Overwrite out_conn with only significant models.
  dbWriteTable(con, 'extra', model_summaries, overwrite = TRUE)
  dbWriteTable(con, 'weights', weights, overwrite = TRUE)
  
  
  dbExecute(con, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbExecute(con, "CREATE INDEX weights_gene ON weights (gene)")
  dbExecute(con, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
  dbExecute(con, "CREATE INDEX gene_model_summary ON extra (gene)")
  
  
  dbDisconnect(con)
  
  print(paste0("Significant genes: ",length(unique(model_summaries$gene))))
  
  ## Return path to database file
  db_path = file.path(saved_db_file_path, output_drv_name)
  return(db_path)
} 
