##############################


# Set working directory and base directory
# setwd("path/to/project/root")


library(purrr);library(limma);library(dplyr);library(tidyverse)




regions = c("dlpfc","dACC","sACC","caudate","hippo","amygdala")


walk(regions, function(region){
  ## Performance and Predictions in testing data
  cis_perf = get(load(paste0("./predictions/CIS/CIS_",region,"_performance_your_testing_genotype_name.RData")))
  cis_perf$Model = "CIS"
  
  if(region != "amygdala"){
    epi_perf = get(load(paste0("./predictions/EpiXcan/EpiXcan_",region,"_performance_your_testing_genotype_name.RData")))
    epi_perf$Model = "EpiXcan"
    
    ## Get common genes and take model based on max adj rsq
    epi.cis_perf       = inner_join(epi_perf, cis_perf, by = c("gene")) %>% 
      dplyr::rename(r_square_adjusted_epi = "r_square_adjusted.x", r_square_adjusted_cis = "r_square_adjusted.y") %>%
      mutate(win_model = case_when(r_square_adjusted_epi > r_square_adjusted_cis ~ "EpiXcan",
                                   r_square_adjusted_cis > r_square_adjusted_epi ~ "CIS",
                                   r_square_adjusted_epi == r_square_adjusted_cis ~ "Same"))
    
    print(table(epi.cis_perf$win_model))
    
    ## Subset performance dataframes with corresponding winning genes
    cis_perf  = cis_perf %>% filter(!gene %in% epi.cis_perf$gene[epi.cis_perf$win_model=="EpiXcan"])
    ## If same keep EpiXcan
    cis_perf  = cis_perf %>% filter(!gene %in% epi.cis_perf$gene[epi.cis_perf$win_model=="Same"])
    
    epi_perf  = epi_perf %>% filter(!gene %in% epi.cis_perf$gene[epi.cis_perf$win_model=="CIS"])
    
    ## save 
    save(cis_perf, file = (paste0("./predictions/CIS/CIS_",region,"your_testing_genotype_name_performance_selected.RData")))
    save(epi_perf, file = (paste0("./predictions/EpiXcan/EpiXcan_",region,"your_testing_genotype_name_performance_selected.RData")))
    
  }else{
    save(cis_perf, file = (paste0("./predictions/CIS/CIS_",region,"your_testing_genotype_name_performance_selected.RData")))
  }
 
  
})





