#!/bin/sh
#$ -N var_cattle
#$ -cwd
#$ -l h_rt=00:25:00
#$ -P roslin_prendergast_cores
#$ -l rl9=true
#$ -l h_vmem=16G

# Run with qsub run_get_variants.cattle.sh <cactus vcf>
# Should run in a few minutes

# Initialise the environment modules
. /etc/profile.d/modules.sh

# Load conda environment and virtualenv
module load anaconda
conda activate tools

VCF=$1
SCRIPTDIR=$2	# directory with run scripts in

${SCRIPTDIR}/get_variants.sh -s 2 -n "BOAEC_D1,Bi_TaA" -v ${VCF} -o var_cattle -k chr_list.txt -p ${SCRIPTDIR}
