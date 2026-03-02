#!/bin/bash


# Get base_name from command line argument
base_name="$1"

# Specify the directory
DIRECTORY="./STEP4.LDproxy_Testing_Geno/LD_output"

output_file="$DIRECTORY/$base_name.vcor"

# Check if the output file already exists
if [[ -f "$output_file" ]]; then
    echo "Skipping $base_name: $output_file already exists."
else
    echo "Processing $base_name..."
    ## Run PLINK2.0
    ./tools/plink2 --pfile ./dataset/g1000/GRCh38/g1000_eur  --ld-snp-list "$DIRECTORY/$base_name".txt --ld-window 1000 --ld-window-kb 500 --ld-window-r2 0.8 --out "$DIRECTORY/$base_name" --r-unphased
fi
