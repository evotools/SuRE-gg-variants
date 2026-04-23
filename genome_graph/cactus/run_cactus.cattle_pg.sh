#!/bin/sh
#$ -N cactus_cow
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=16G
#$ -P roslin_prendergast_cores
#$ -l rl9=true

# Run with
# qsub -pe sharedmem 16 -R y run_cactus.cattle_pg.sh

# Need chr_list in the directory

# Initialise the environment modules
. /etc/profile.d/modules.sh

SCRIPTDIR=$1

# Load virtualenv for cactus
source /exports/cmvm/eddie/eb/groups/prendergast_roslin/Lindsey/programs/cactus/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate

# Run cactus specifying the chromosomes including X and Y and keeping other contigs
${SCRIPTDIR}/cactus_pangenome.sh -s seqfile.cattle.txt -o cattle-Tann-pg -r ARS_UCD2_Tann -a "ARS_UCD2_Tann Bi_83_A BOAEC_D1" -c "--refContigs 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 X Y CR940347.1 CR940349.1 CR940350.1 CR940351.1 CR940348.1 CR940352.2 CR940353.1 CR940346.1 --otherContig chrOther"
 
${SCRIPTDIR}/extract_walks.sh -n "Bi_83_A,BOAEC_D1" -k chr_list.txt -r ARS_UCD2_Tann -o cattle-Tann-pg
