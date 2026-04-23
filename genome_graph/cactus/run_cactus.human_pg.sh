#!/bin/sh
#$ -N cactus_hum
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=16G
#$ -P roslin_prendergast_cores
#$ -l rl9=true

# Run with
# qsub -pe sharedmem 16 -R y run_cactus.human_pg.sh

# Initialise the environment modules
. /etc/profile.d/modules.sh

SCRIPTDIR=$1

# Load virtualenv for cactus
source /exports/cmvm/eddie/eb/groups/prendergast_roslin/Lindsey/programs/cactus/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate

${SCRIPTDIR}/cactus_pangenome.sh -s seqfile.hum.txt -o human-pg -r CHM13 -a "CHM13 NA12878" -c "--refContigs chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 --otherContig chrOther"

${SCRIPTDIR}/extract_walks.sh -n "NA12878" -k chrlist.hum.txt -r CHM13 -o human-pg