#!/bin/sh
#$ -N var_human
#$ -cwd
#$ -l h_rt=01:00:00
#$ -P roslin_prendergast_cores
#$ -l rl9=true
#$ -l h_vmem=16G

# Run with qsub run_get_variants.human.sh <cactus vcf>
# Should have list of IDs of chromosomes to keep in chr_list.txt in wd 
# Should run in a few minutes

# Initialise the environment modules
. /etc/profile.d/modules.sh

# Load conda environment and virtualenv
module load anaconda
conda activate tools

VCF=$1
SCRIPTDIR=$2	# directory with run scripts in

${SCRIPTDIR}/get_variants.sh -s 1 -n "NA12878" -v ${VCF} -o var_NA12878 -k chr_list.txt -p ${SCRIPTDIR}