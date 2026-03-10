# coTWAS - co-expression TWAS pipeline

## Overview

This pipeline trains and applies multiple gene expression prediction models for brain tissues using RNA-seq and genotype data. It combines cis-regulatory (CIS, EpiXcan) and trans-regulatory (MODULE, INGENE) approaches to predict gene expression, then performs co-expression-based Transcriptome-Wide Association Studies (coTWAS).


## Citation

If you use this pipeline, please cite the corresponding manuscript:

Rossi et al., [INSERT DOI]

This pipeline integrates methods from PrediXcan, EpiXcan, MetaXcan, MAGMA, GCTA-COJO, and elastic-net modeling (glmnet). Appropriate citation of these tools is required in derivative work.

R packages used and version information are provided in the Key Resources Table in the Supplementary Materials of the associated manuscript and in the 'README_Packages_and_Versions' file within this repository.


## Main Objective

The pipeline builds predictive models that estimate gene expression levels from genotype data across multiple brain regions (DLPFC, Caudate, Hippocampus, Amygdala, dACC, sACC). It uses four complementary approaches:

- **CIS**: cis-eQTL models using variants within 1Mb of genes
- **EpiXcan**: Epigenetic annotation-informed cis-eQTL models 
- **MODULE**: Co-expression network-based model using ME of coexpressed genes
- **INGENE**: Co-expression network-based model using imputed expression from co-expression partners

## Pipeline Structure

The pipeline consists of 13 sequential steps:

### STEP1: Process Data
- Preprocesses RNA-seq data (SummarizedExperiment format)
- Removes outlier samples and genes
- Calculates cell type proportions
- **Input**: Raw RSE data, genomic eigenvariables
- **Output**: Cleaned expression data (`outlierRemoved.rds`)

### STEP2: Training CIS
- Performs cis-eQTL analysis
- Trains elastic net models for each gene using variants in cis window
- Creates SNP annotation lookup for EpiXcan
- **Input**: Expression data, genotype data by chromosome
- **Output**: Model weights and summaries per chromosome

### STEP3: Training MODULE
- Maps genes to co-expression network modules
- Identifies co-expression eQTLs (coeQTLs)
- Trains network-based prediction models
- **Input**: Expression data, genotype data, network modules
- **Output**: Network-specific model weights and summaries

### STEP4: LDproxy Testing Geno
- Tests genotype matching using LD proxies
- Validates genotype compatibility between training and testing datasets

### STEP5: Match Geno
- Matches genotypes between training and testing datasets
- Prepares genotype data for prediction

### STEP6: Apply EpiXcan CIS MODULE Models
- Applies trained CIS, EpiXcan, and MODULE models to testing data
- Generates predicted expression values
- Uses MetaXcan framework 

### STEP7: Training INGENE
- Trains INGENE models using imputed expression from CIS/EpiXcan
- Uses co-expression network partners as features
- **Input**: Imputed expression from STEP6, network modules
- **Output**: INGENE model weights and summaries

### STEP8: Validate CIS Predictions Testing Data
- Validates prediction performance on testing data
- Calculates correlation and R² metrics

### STEP9: Apply INGENE
- Applies trained INGENE models to testing data
- Generates INGENE-based expression predictions

### STEP10: Average Network Predictions
- Averages INGENE and MODULE predictions across multiple networks for each gene
- Combines predictions from different co-expression networks
- **Output**: Averaged network predictions

### Cross-Training Selection: Replicable Genes (Optional)
- **Purpose**: Identifies replicable trans-regulatory genes (MODULE and INGENE) across multiple training datasets
- **When to use**: If you have trained models on more than one dataset (e.g., GTEx, LIBD, CMC), this step identifies genes with consistent predictions across datasets
- **Method**: Compares predictions from different training datasets on the same testing data and calculates gene-specific correlations
- **Input**: Averaged network predictions from multiple training datasets (e.g., GTEx-trained vs LIBD-trained predictions)
- **Output**: Lists of replicable genes per brain region (correlation > threshold) saved as RDS files
- **Usage**: Run `Cross_Training_Selection_Replicable_Genes/MODULE_selection.R` and `Cross_Training_Selection_Replicable_Genes/INGENE_selection.R` after completing STEP10 for all training datasets

### STEP11: Combine CIS Trans Predictions
- Combines predictions from CIS/EpiXcan and MODULE/INGENE methods
- Fits linear models to optimize combination weights
- **Output**: Combined prediction models

### STEP12: MLE
- Maximum Likelihood Estimation analysis

### STEP13: coTWAS Analysis
- Performs co-expression-based Transcriptome-Wide Association Studies
- Combines predictions across methods for association testing
- **Input**: Combined predictions
- **Output**: coTWAS association results

## Prerequisites

- **R** (version 4.0 or higher)
- **Required R packages**: See individual step scripts for package requirements and in the 'README_Packages_and_Versions' in this folder
- **Python** (version 3.0 or higher) (for MetaXcan scripts in STEP6)
- **Perl** (for MetaXcan wrapper scripts)

###  Run Pipeline Sequentially

Execute steps in order:

```r
# STEP1: Preprocess data
source("STEP1.Process_Data/1.RemoveOutliers.R")
source("STEP1.Process_Data/2.cleanCovs.R")

# STEP2: Train CIS models
source("STEP2.Training_cis/1.subsetGenotype.R")
source("STEP2.Training_cis/2.cisEQTL.R")
source("STEP2.Training_cis/3.Training.R")
source("STEP2.Training_cis/3.Training_EpiXcan.R")

# STEP3: Train MODULE models
source("STEP3.Training_MODULE/1.prepare_mapping.R")
source("STEP3.Training_MODULE/2.mapping.R")
source("STEP3.Training_MODULE/3.coeQTL.R")
source("STEP3.Training_MODULE/4.Training.R")

# STEP10: Average network predictions (after applying models)
# ... (apply models and average predictions)

# Optional: Cross-Training Selection (if using multiple training datasets)
# source("Cross_Training_Selection_Replicable_Genes/MODULE_selection.R")
# source("Cross_Training_Selection_Replicable_Genes/INGENE_selection.R")

# Continue with remaining steps...
```

### 4. Key Configuration Parameters

Common parameters to adjust across steps:

- **Brain regions**: `c("dlpfc", "caudate", "hippo", "amygdala", "dACC", "sACC")`
- **Chromosomes**: `1:22`
- **Cis window**: `1e6` (1Mb)
- **Cross-validation folds**: `4`
- **Elastic net alpha**: `0.5` 
- **Number of cores**: Adjust based on your system

## Output Structure

```
output/
├── STEP2.Training_cis/output/
│   ├── summary/          # Model performance summaries
│   └── weights/           # SNP weights
├── STEP3.Training_MODULE/output/
│   ├── summary/          # Network model summaries
│   └── weights/          # Network model weights
├── STEP7.Training_INGENE/output/
│   ├── summary/          # INGENE model summaries
│   └── weights/          # INGENE model weights
├── predictions/          # Predicted expression values
└── coTWAS/              # Final association results
```

## Notes

- Each step checks for existing output files and can skip completed analyses
- The pipeline supports parallel processing (adjust `n_cores` in each step)
- Some steps require significant computational resources and memory
- The pipeline uses elastic net regression (alpha=0.5) for all model training


