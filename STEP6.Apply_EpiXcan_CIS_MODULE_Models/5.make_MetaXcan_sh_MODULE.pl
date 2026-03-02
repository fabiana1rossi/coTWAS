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
my $output_directory = "./STEP6.Apply_EpiXcan_CIS_MODULE_Models/MetaXcan_sh_scripts";
my $metaxcan_path    = "../../tools/Metaxcan/software";
my $main_path        = "../../PredictDB/module/";
my $data_path        = "./genotype/MODULE.matched.geno";
my $results_path     = "./predictions/";

# List of genotype directories
my @genotype_files = ("your_testing_genotype_prefix"); 
my @targets        = ("your_testing_genotype_name");
my @regions        = ("dlpfc");

# List of networks
my @networks = ("caudate", "caudate__.1.25", "caudate__25.50", "caudate__50.100", "dentate",
                "dentate.noQSVAremoved", "dlpfc", "dlpfc__.1.6", "dlpfc__6.25", "dlpfc__25.50", 
                "dlpfc__50.100", "Fromer2016_control", "Gandal2018", "Gandal2018PE", "Gandal2018PE_cs",
                "HartlALL", "HartlBGA", "HartlBRNACC", "HartlBRNAMY", "HartlBRNCBH", "HartlBRNCBL",
                "HartlBRNCDT", "HartlBRNCTX", "HartlBRNCTXB24", "HartlBRNCTXBA9", "HartlBRNHIP",
                "HartlBRNHYP", "HartlBRNPUT", "HartlBRNSNA", "HartlBROD", "HartlCEREBELLUM",
                "HartlCTX", "HartlNS.SCTX", "HartlSTR", "HartlSUBCTX", "HartlWHOLE_BRAIN", 
                "hippo", "hippo.QSVAremoved", "hippo__.1.6", "hippo__6.25", "hippo__25.50", 
                "hippo__50.100", "Li2018", "Pergola2019", "Pergola2020", "Radulescu2020",
                "Walker2019", "Werling2020");


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
    
            # Ensure the RESULTS directory structure exists
            my $results_dir = "$results_path/${target}/MODULE";
            unless (-d $results_dir) {
                make_path($results_dir) or die "Failed to create directory $results_dir: $!";
                print "Directory created: $results_dir\n";
            }

        foreach my $network (@networks) {
            # Loop over genotype identifiers
            foreach my $genotype_file (@genotype_files) {

                # Declare the variable $geno_id (this is not necessary, but if you need it later, use it)
                my $geno_id = "$genotype_file";
                # You already have $network and $region from the outer loops
                
                # Generate the script path
                my $script_path = "${output_directory}/MODULE_${region}_MetaXcan_script_${target}_${network}.sh";
                print "    Script path: $script_path\n";

                # Open the output file for the bash script
                open my $script_fh, '>', $script_path or die "Unable to open $script_path: $!";
                
                # Write the bash script content
                print $script_fh <<"EOF";
#!/bin/bash

export METAXCAN=$metaxcan_path
export MAIN=$main_path
export DATA=$data_path
export RESULTS=$results_path

printf "Predict expression: $target \\n\\n"
python3 \$METAXCAN/Predict.py \\
--model_db_path \$MAIN/${region}/${region}__LIBD.chr1_22.MODULE_EUR_${network}.db \\
--model_db_snp_key rsid \\
--vcf_genotypes \$DATA/${genotype_file}_${region}_MODULE.matched.vcf \\
--vcf_mode genotyped \\
--force_colon \\
--prediction_output \$RESULTS/MODULE/MODULE_${region}_${network}_predicted.txt.gz \\
--prediction_summary_output \$RESULTS/MODULE/MODULE_${region}_${network}_summary.txt.gz \\
--verbosity 9 \\
--throw

EOF

                close $script_fh;

                # Set the permissions to make the script executable
                chmod 0755, $script_path;
            }
        }
    }
}

print "Scripts have been generated in '$output_directory'.\n";
