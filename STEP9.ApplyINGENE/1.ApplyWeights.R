#####################################################################
# INGENE GTEx Prediction Pipeline
#####################################################################
# Applies trained INGENE models to GTEx data to generate predictions
#####################################################################

## Set working directory
# Set working directory to project root
# setwd("path/to/project/root")

base_dir = getwd()


source('./STEP9.ApplyINGENE/apply_weights_helper_functions.R')

# Required packages
required_packages = c(
    "purrr", "parallel", "dplyr", "limma", 
    "data.table"
)

# Install and load required packages
for(pkg in required_packages) {
    if(!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
        library(pkg, character.only = TRUE)
    }
}

#################################################
## Configuration
#################################################
testing_genotype        = ""
testing_genotype_folder = ""

CONFIG = list(
    input = list(
        # Data types
        cis_types = c("EpiXcan","CIS"),
        regions = c("dlpfc"),
        
        # Input files
        weights_dir   = "./STEP7.Training_INGENE/output/weights/",
        summaries_dir = "./STEP7.Training_INGENE/output/summary/",
        
        # CIS/EpiXcan predictions
        cis     = "./predictions/%s/CIS/CIS_%s_%s_predicted.txt.gz",
        epixcan = "./predictions/%s/EpiXcan/EpiXcan_%s_%s_predicted.txt.gz"
        
        # CIS/EpiXcan training genes
       # training_data = "./predictions/LIBD_training/CIS_Epi_Imputed.RData"
    ),
    
    output = list(
        base_dir = sprintf('./predictions/%s/INGENE/',testing_genotype_folder),
        file_prefix = "INGENE_"
    ),
    
    compute = list(
        n_cores = 20
    )
)

# Create output directories
for(region in CONFIG$input$regions) {
    for(cis in CONFIG$input$cis_types) {
        dir.create(file.path(CONFIG$output$base_dir, cis, region), 
                  recursive = TRUE, showWarnings = FALSE)
    }
}

#################################################
## Data Loading Functions
#################################################

#' Load imputed expression data
#' @param config Configuration list
#' @param data_type Type of data to load ("epixcan" or "cis")
#' @return List of imputed expression data by brain region
load_imputed_data = function(config, region, data_type) {
  
    if(region == "amygdala" && data_type == "EpiXcan"){return(as.data.frame(matrix(NA,ncol = 1, nrow = 1)))} ## No EpiXcan annotations/training for this region
  
    input_config = sprintf(config$input[[tolower(data_type)]],testing_genotype_folder, region, testing_genotype)
    
    # Load or create predictions
    predictions = fread(input_config) %>% as.data.frame()
    
    rownames(predictions) = predictions$IID
   return(predictions)
}

#################################################
## Main Pipeline
#################################################

main <- function() {
    start_time = Sys.time()
    
    # Process each CIS type
    for(cis_type in CONFIG$input$cis_types) {
        message("Processing ", cis_type)
        if(cis_type=="EpiXcan"){regions=c("dlpfc")}else{regions="dlpfc"}
      
        # Process each region
        for(region in CONFIG$input$regions) {  
            message("Processing region: ", region)
            
            # Load predictions
            predictions = load_imputed_data(CONFIG, region, ifelse(cis_type == "CIS", "CIS", "EpiXcan"))
            if(all(dim(predictions) == c(1,1))) {
              warning("No predictions found for ", cis_type)
              next
            }
            
          
            # # Load training winning genes
            # selected_genes = load(training_data)
            # selected_genes = selected_genes[[tolower(cis_type)]]
            
            
            # Load summaries
            summaries_file_path = file.path(CONFIG$input$summaries_dir, region)
            if(!file.exists(summaries_file_path)) {
              warning("Missing summaries folder: ", summaries_file_path)
              next
            }
            
            
            # Load weights
            weights_file_path = file.path(CONFIG$input$weights_dir, region)
            if(!file.exists(weights_file_path)) {
                warning("Missing weights folder: ", weights_file_path)
                next
            }
            
            ## Get network names
            networks = strsplit2(list.files(weights_file_path),"_weights")[,1]
            
            if(is.null(networks)) next
            
            # Get network predictions
            walk(unique(networks), function(net) {
                message(sprintf("Processing network %d/%d", which(unique(networks)==net), length(unique(networks))))
              
                ## open weight file
                net_weights   = read.delim(file.path(weights_file_path,paste0(net,"_weights.txt")))
                net_summaries = read.delim(file.path(summaries_file_path,paste0(net,"_summary.txt")))
                
                ## Filter only significant genes
                net_summaries  = net_summaries[which(net_summaries$adj_rsq_gene >= 0.01 & pval_gene < 0.05),]
                net_weights    = net_weights[which(net_weights$gene %in% net_summaries$gene_id),]
                
                if(nrow(net_weights) == 0) return(NULL)
                
                output_path = file.path(CONFIG$output$base_dir, cis_type, region,
                                      paste0(CONFIG$output$file_prefix,region, "_",  
                                            net, "_predictions.txt"))
                
                computePredFromWeights(
                    test.pred   = predictions,
                    mod.weights = net_weights,
                    network     = net,
                    output_file = output_path,
                    n_cores     = CONFIG$compute$n_cores
                )
            })
          }
    }
    
    end_time = Sys.time()
    message("Pipeline completed in ", difftime(end_time, start_time, units="mins"), " minutes")
}

# Run pipeline
main()
