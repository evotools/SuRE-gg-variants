#!/bin/sh
#$ -N fragpath
#$ -cwd
#$ -l h_vmem=24G
#$ -l h_rt=24:00:00
#$ -P roslin_prendergast_cores
#$ -pe sharedmem 2
#$ -t 1-29
#$ -R y

BTPATH=$1	# path to BOAEC-D1 walks directory
BIPATH=$2	# path to Bi-83-A walks directory
BCPATH=$3
BCFILE=$4	# filename suffix within bcpath directory

bclist="${SGE_TASK_ID}.barcodes.txt"
if ! [ -f "${bclist}" ]; then
	echo "No barcodes for chr ${SGE_TASK_ID}"
	exit 0
fi

allbc="${BCPATH}/${SGE_TASK_ID}.${BCFILE}"
echo "Using barcode file ${allbc}"
if ! [ -f "${allbc}" ]; then
	echo "Expected barcodes file ${allbc} doesn't exist"
	exit 1
fi

# Get the barcodes and split into libraries
gawk 'NR==FNR{a[$1]; next} $1 in a' ${bclist}  <(zcat ${allbc}) > ${SGE_TASK_ID}.bc_to_check.txt
grep 'Bi-83-A' ${SGE_TASK_ID}.bc_to_check.txt | sort -k1,1 > Bi.${SGE_TASK_ID}.bc_to_check.txt
grep 'BOAEC-D1' ${SGE_TASK_ID}.bc_to_check.txt | sort -k1,1 > Bt.${SGE_TASK_ID}.bc_to_check.txt

btwalks="${BTPATH}/${SGE_TASK_ID}.BOAEC_D1.walks.txt"
biwalks="${BIPATH}/${SGE_TASK_ID}.Bi_83_A.walks.txt"

# Run the script
./find_frag_path.pl ${btwalks} Bt.${SGE_TASK_ID}.bc_to_check.txt Bt.${SGE_TASK_ID}.frag_paths.txt
./find_frag_path.pl ${biwalks} Bi.${SGE_TASK_ID}.bc_to_check.txt Bi.${SGE_TASK_ID}.frag_paths.txt
