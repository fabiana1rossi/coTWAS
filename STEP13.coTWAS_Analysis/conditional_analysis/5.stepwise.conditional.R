# Load necessary libraries
library(coco)

# Set base directory for conditional analysis results
# base_dir = "path/to/conditional/results/"
base_dir = "./result/conditional.results/"

##### load Loci
load(file.path(base_dir, "Loci.RData"))
to.keep = names(table(Loci$Region))[table(Loci$Region) > 1]
check = gsub("\\.RData","",
             limma::strsplit2(grep("results_",list.files(base_dir),value = T),"_")[,2])
to.keep = to.keep[!to.keep %in% check]

purrr::walk(to.keep,function(locus){
  print(locus)
  # Step 1: Load gene-level summary statistics (significant genes only, pFDR <= 0.01)
  load(paste0(base_dir, "summary/summary.input_",locus,".RData"))
  
  # Step 2: Load custom LD matrix (correlation matrix among gene predictors)
  load(paste0(base_dir, "corr.matrix/cor.matrix_",locus,".RData"))
  
  # remove genes with all NA in the LD matrix
  to.remove = apply(avg_correlation_matrix,2,function(x)all(is.na(x)))
  to.remove = names(to.remove)[to.remove]
  
  # Ensure that the LD matrix and gene_stats are aligned
  common_genes <- intersect(rownames(avg_correlation_matrix), summary$ID)
  common_genes = common_genes[!common_genes %in% to.remove]
  ld_matrix <- avg_correlation_matrix[common_genes, common_genes]
  gene_stats <- summary[match(common_genes, summary$ID), ]
  colnames(gene_stats)[1] = "rsid"
  gene_stats$N = as.numeric(gene_stats$N)
  gene_stats = gene_stats[!duplicated(gene_stats$beta),]
  
  #### order based on pval
  gene_stats = gene_stats[order(gene_stats$p),]
  ld_matrix <- ld_matrix[gene_stats$rsid, gene_stats$rsid]
  
  # Step 3: Format input for CoCo
  vars = gene_stats$GVAR * (gene_stats$N) * gene_stats$SE^2 * (gene_stats$N - 1) + gene_stats$GVAR * (gene_stats$N) * gene_stats$beta^2
  var_y = (median(vars/(gene_stats$N - 1)))
  coco_input <- prep_dataset_coco(data_set = gene_stats,ld_matrix = ld_matrix,var_y = var_y)
  coco_input$data_set$var = coco_input$data_set$var[[1]]
  coco_input$data_set$neff = coco_input$data_set$neff[[1]]
  coco_input$hwe_diag = coco_input$hwe_diag[[1]]
  coco_input$hwe_diag_outside = coco_input$hwe_diag_outside[[1]]
  
  # Step 4: Run forward stepwise conditional analysis
  # stepwise_joint_results <- stepwise_coco(coco_input,p_value_threshold = 0.01,joint = T,return_all_betas = F,max_iter = 2184)
  stepwise_cond_results <- stepwise_coco(coco_input,p_value_threshold = 1,joint = F,return_all_betas = F,max_iter = nrow(gene_stats))
  
  # Step 5: Adjust CoCo betas to unit variance using GVAR
  # gene_stats_joint = merge(gene_stats,stepwise_joint_results$stepwise_summary,by = "rsid")
  # gene_stats_joint$Adjusted_Beta <- gene_stats_joint$beta_new * sqrt(gene_stats_joint$GVAR)
  # gene_stats_joint$Adjusted_OR <- exp(gene_stats_joint$Adjusted_Beta)
  gene_stats_cond = merge(gene_stats,stepwise_cond_results$stepwise_summary,by = "rsid")
  # gene_stats_cond$p.fdr = p.adjust(gene_stats_cond$p_new, method = "BH")
  gene_stats_cond$Adjusted_Beta <- gene_stats_cond$beta_new * sqrt(gene_stats_cond$GVAR)
  gene_stats_cond$Adjusted_OR <- exp(gene_stats_cond$Adjusted_Beta)
  
  # Step 6: Save results
  save(gene_stats_cond, 
       file = paste0(base_dir, "results_",locus,".RData"))
  
})



