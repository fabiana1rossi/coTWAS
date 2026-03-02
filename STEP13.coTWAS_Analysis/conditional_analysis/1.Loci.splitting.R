
# Load gene-level association results
load("PGC3_snellius/coTWAS/result/output_all.RData")
result$tissue.gene = paste0(result$gene,".",result$tissue)
result$chr = gsub("chr","",as.character(result$seqnames))
result.noMHC = result[result$MHC %in% "noMHC",]

# Ensure correct ordering
genes <- result.noMHC[order(result.noMHC$chr, result.noMHC$start), ]

# Set maximum distance threshold (e.g., 1 Mb)
distance_threshold <- 1e6

# Initialize region assignment
genes$Region <- NA
current_region = 0

for (chr in c(as.character(1:22),"X","Y")){
  chr_genes <- genes[genes$chr == chr, ]
  current_region = current_region + 1
  chr_genes$Region <- current_region
  
  for (i in 2:nrow(chr_genes)) {
    distance <- chr_genes$start[i] - chr_genes$end[i - 1]
    if (distance > distance_threshold) {
      current_region <- current_region + 1
    }
    chr_genes$Region[i] <- current_region
  }
  
  genes[genes$chr == chr, ] <- chr_genes
  
}
genes$Region = as.character(genes$Region)

# set MHC region
genes.MHC = result[result$MHC %in% "MHC",]
genes.MHC$Region = "MHC"

# Save region-annotated file
Loci = rbind(genes,genes.MHC)
save(Loci, file = "PGC3_snellius/coTWAS/result/conditional.results/Loci.RData")


