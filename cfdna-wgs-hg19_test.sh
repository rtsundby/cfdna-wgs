#! /bin/bash

module load snakemake
#module load singularity
cd /data/ShernLiquidBx/cfdna-wgs

# Necessary to run conda snakemake command in shell script
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

##################################################################
###                          TESTING                           ###
##################################################################
## dry-run 
# snakemake -s workflow/int_test.smk --configfile config/int_test_biowulf.yaml --profile profile/ -np
## actual run
# snakemake -s workflow/int_test.smk --configfile config/int_test_biowulf.yaml --profile profile/

##################################################################
###                        MPNST dataset                       ###
##################################################################
## dry-run 
# snakemake -s workflow/hg19-slim.smk --configfile config/nf1_hg19_test.yaml --profile profile/ -np

## actual run
snakemake -s workflow/hg19-slim.smk --configfile config/nf1_hg19_test.yaml --profile profile/

## START SCRIPT:
# sbatch --mem=1g --cpus-per-task=2 --time=72:00:00 --mail-type=BEGIN,TIME_LIMIT_90,END /data/ShernLiquidBx/cfdna-wgs/cfdna-wgs-hg19.sh
