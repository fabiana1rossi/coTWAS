library(data.table)
library(dplyr)
library(parallel)
library(limma)
library(purrr)

# Set working directory to project root
# setwd("path/to/project/root")



### Params
target_file_path = c("./genotype/")
site               = "your_geno_bim" # LIBD_TopMed.geno.maf.EUR.chr1_22.noindels.overlap.GTEX.bim
name_site          = "your_geno_name"

## Extract target genotype bim SNP file
target           =  fread(file.path(target_file_path,site))


## g1000 bim file
bim_g1000       = read.delim("./dataset/g1000/GRCh38/g1000_eur.pvar") ## from 1000G




print("start")
if(!file.exists(paste0("./STEP4.LDproxy_Testing_Geno/LD_output/MODULE.", name_site, "_LDtable.RData"))){
  
  
  # Read in output
  LD0 = fread(paste0("./STEP4.LDproxy_Testing_Geno/LD_output/",name_site,"_mismatch.vcor"))
  LD = LD0[LD0$ID_A != LD0$ID_B,]
  
  # Use mclapply to parallelize within the task, depending on the number of cores
  LD_chunk = mclapply(unique(LD$ID_A), function(x) {
    print(paste0(which(unique(LD$ID_A) == x),"/",length(unique(LD$ID_A))))
    df = LD[LD$ID_A == x, ]
    
    proxy = df$ID_B[df$ID_B %in% target$V2]
    df = df[order(abs(df$UNPHASED_R), decreasing = TRUE), ]
    df = df[df$ID_B %in% proxy,]
    
    if (length(proxy)) {
      df.g1000 = bim_g1000[which(bim_g1000$ID %in% proxy),]
      table.match = inner_join(df,df.g1000, by=c("ID_B"="ID"))
      table.match = table.match %>% dplyr::select(ID_A,MAJ_A,ID_B,MAJ_B,UNPHASED_R,REF,ALT) %>%
        dplyr::rename(ref = "ID_A",target.LD="ID_B",REF.g1000 = "REF",ALT.g1000="ALT")
      
    } else {
      table.match = data.frame(ref = x, MAJ_A=NA, target.LD = x, MAJ_B=NA, UNPHASED_R=NA,REF.g1000 = NA,ALT.g1000=NA)
    }
    
    return(table.match)
  }, mc.cores = 8)  # Adjust this based on the number of cores available 
  
  # Combine the results from the chunk
  LD = do.call('rbind', LD_chunk)
  LD = LD %>% mutate(check = case_when(ref == target.LD ~ 'same', ref != target.LD ~ 'different'))
  
  # Save the result for this task
  save(LD, file = paste0("./STEP4.LDproxy_Testing_Geno/LD_output/MODULE.", name_site, "_LDtable.RData"))
  
  # Print completion message
  cat("Completed task", name_site, "\n")
}else{ # Print completion message
  cat("ALREDY DONE task", name_site, "\n")
}

