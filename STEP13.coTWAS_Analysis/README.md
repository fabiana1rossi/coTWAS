# coTWAS Analysis Pipeline

This directory contains scripts for performing coTWAS (co-expression-based Transcriptome-Wide Association Studies) analysis. The pipeline consists of three main steps:

1. Preprocessing
2. Analysis
3. Results Processing

## Prerequisites

- R (version 4.0 or higher)
- Required R packages:
  - furrr
  - data.table
  - purrr
  - rcompanion
  - broom
  - dplyr
  - mgsub
  - metafor
  - metap

## Directory Structure

```
STEP12.coTWAS/
├── config.R           # Configuration file with paths and parameters
├── utils.R           # Utility functions
├── 1.coTWAS_preprocessing.R
├── 2.PGC_TWAS_Analyses.R
├── 3.coTWAS_output.processing.R
└── README.md
```

## Usage

### 1. Preprocessing

The preprocessing step prepares the input data for analysis:

```bash
Rscript 1.coTWAS_preprocessing.R <site_index>
```

Where:
- `site_index`: Index of the site to analyze (numeric)

### 2. Analysis

The analysis step performs the coTWAS analysis:

```bash
Rscript 2.PGC_TWAS_Analyses.R <tissue> <site_index>
```

Where:
- `tissue`: Tissue to analyze (e.g., "dlpfc", "hippo", etc.)
- `site_index`: Index of the site to analyze (numeric)

### 3. Results Processing

The processing step combines and analyzes results across sites:

```bash
Rscript 3.coTWAS_output.processing.R <tissue>
```

Where:
- `tissue`: Tissue to process (e.g., "dlpfc", "hippo", etc.)

## Configuration

Before running the scripts, make sure to:

1. Update the paths in `config.R` to match your system
2. Ensure all required input files are present
3. Create necessary output directories

## Input Files Required

- Metadata file (specified in `config.R`)
- Common genes file (specified in `config.R`)
- Unique genes file (specified in `config.R`)
- Gene information file (specified in `config.R`)
- Prediction files in the specified directories

## Output Files

The pipeline generates the following output files:

1. Preprocessing:
   - `{site}_{tissue}_predictions.common.RData`
   - `{site}_{tissue}_predictions.unique.RData`

2. Analysis:
   - `{site}_{tissue}_coTWAS.Analysis.RData`

3. Processing:
   - `output_{tissue}.RData`

## Error Handling

The scripts include error handling for:
- Missing input files
- Invalid command line arguments
- Data validation
- File I/O operations

## Notes

- The pipeline uses parallel processing where appropriate
- Memory usage is optimized for large datasets
- Progress messages are displayed during execution 