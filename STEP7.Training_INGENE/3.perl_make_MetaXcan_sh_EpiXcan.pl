#!/usr/bin/perl

# This script automates the generation of shell scripts for running MetaXcan predictions.
# It ensures that necessary directories exist before execution and creates output scripts
# that define environment variables and execute the MetaXcan prediction tool.
#
# The script processes a list of genotype targets and regions, generating corresponding
# shell scripts for each combination. It also ensures the results directory structure exists
# before writing the scripts.

use strict;
use warnings;
use File::Path qw(make_path);

# Define variables
my $output_directory = "./STEP7.Training_INGENE/MetaXcan_sh_scripts";
my $metaxcan_path    = "../../tools/Metaxcan/software";
my $main_path        = "./STEP2.Training_CIS/output_GTEx_EpiXcan/database/";
my $data_path        = "./genotype/EpiXcan.matched.geno";
my $results_path     = "./predictions";
my $training_dataset_folder = "GTeX.v9";
my $training_dataset = "GTEx";
my $outdir           = "GTEx_training";
my @genotype_files = ("GTeX.v9.phased.geno.maf.bfile.updated.EUR"); # ("training_genotype")
my @targets        = ("GTeX.v9.phased.geno.maf.bfile.updated.EUR"); # ("training_genotype")
my @regions        = ("acc","dlpfc","caudate","hippo");

# Create the output directory if it doesn't exist
unless (-d $output_directory) {
    mkdir $output_directory or die "Failed to create directory $output_directory: $!";
    print "Directory created: $output_directory\n";
} else {
    print "Directory already exists: $output_directory\n";
}

# Loop over genotype identifiers
foreach my $target (@targets) {
    print "Processing target: $target\n";
    foreach my $region (@regions) {
        print "  Processing region: $region\n";
        foreach my $genotype_file (@genotype_files) {
            print "    Processing genotype file: $genotype_file\n";
            
            # Ensure the RESULTS directory structure exists
            my $results_dir = "$results_path/${outdir}/${training_dataset_folder}/EpiXcan";
            unless (-d $results_dir) {
                make_path($results_dir) or die "Failed to create directory $results_dir: $!";
                print "Directory created: $results_dir\n";
            }
            
            # Construct the full path for the sh script
            my $script_path = "${output_directory}/EpiXcan_${region}_MetaXcan_script_${target}.sh";
            print "    Script path: $script_path\n";

            # Open the output file for the sh script
            open my $script_fh, '>', $script_path or die "Unable to open $script_path: $!";
            print $script_fh "#!/bin/bash\n\n";
            print $script_fh "export METAXCAN=$metaxcan_path\n";
            print $script_fh "export MAIN=$main_path\n";
            print $script_fh "export DATA=$data_path\n";
            print $script_fh "export RESULTS=$results_path\n\n";
            print $script_fh "printf \"Predict expression: ${target} in ${region} \\\n\\n\"\n\n";
            print $script_fh "python3 \$METAXCAN/Predict.py \\\n";
            print $script_fh "--model_db_path \$MAIN/${region}/EpiXcan_${region}_${training_dataset}.db \\\n";
            print $script_fh "--model_db_snp_key rsid \\\n";
            print $script_fh "--vcf_genotypes \$DATA/${genotype_file}_${region}_EpiXcan.matched.vcf \\\n";
            print $script_fh "--vcf_mode genotyped \\\n";
            print $script_fh "--force_colon \\\n";
            print $script_fh "--prediction_output \$RESULTS/${outdir}/${training_dataset_folder}/EpiXcan/EpiXcan_${region}_${training_dataset}_predicted.txt.gz \\\n";
            print $script_fh "--prediction_summary_output \$RESULTS/${outdir}/${training_dataset_folder}/EpiXcan/EpiXcan_${region}_${training_dataset}__summary.txt.gz \\\n";
            print $script_fh "--verbosity 9 \\\n";
            print $script_fh "--throw\n\n";
            
            close $script_fh;
            print "    Script generated: $script_path\n";

            # Change permissions
            chmod 0755, $script_path or die "Failed to chmod $script_path: $!";
        }
    }
}

print "All scripts have been generated in '$output_directory'.\n";
