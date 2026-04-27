#!/bin/sh
#$ -N refpos
#$ -cwd
#$ -l h_vmem=24G
#$ -l h_rt=24:00:00
#$ -P roslin_prendergast_cores
#$ -pe sharedmem 2
#$ -t 1-29
#$ -R y

# 1 is path to subgraphs
# 2 is path to barcodes list
# 3 is path to vg


bclist="${SGE_TASK_ID}.barcodes.txt"
if ! [ -f "${bclist}" ]; then
	echo "No barcodes for chr ${SGE_TASK_ID}"
	exit 0
fi

subgraph="${1}/${SGE_TASK_ID}.vg"
echo "Using subgraph $subgraph"
if ! [ -f "${subgraph}" ]; then
	echo "Expected subgraph file ${subgraph} doesn't exist"
	exit 1
fi


allbc="${2}/${SGE_TASK_ID}.all_cattle.multi.signif_frag_positions.txt.gz"
echo "Using barcode file ${allbc}"
if ! [ -f "${allbc}" ]; then
	echo "Expected barcodes file ${allbc} doesn't exist"
	exit 1
fi

VGPATH=$3

./project_to_ref.sh -l ${bclist} -b ${allbc} -g ${subgraph} -v ${VGPATH} -p ${SGE_TASK_ID}
