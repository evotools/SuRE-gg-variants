#!/usr/bin/env bash

# extract_walks.sh

# Extract chromosome subgraphs from the all-refs pangenome
# From these, extract all walks for each subgraph, then get them for each sample and put in sample subdirectory

VGPATH="vg"
OPTSTRING="n:k:r:o:v:"
while getopts ${OPTSTRING} opt; do
	case $opt in
	n)
		NAMES=${OPTARG}		# comma-separated list of sample names
		;;
	k)
		KEEP=${OPTARG}		# list of names of chromosomes to keep
		;;
	r)
		REFID=${OPTARG}
		;;
	o)
		OUTPUT=${OPTARG}	# directory used for cactus_pangenome.sh output
		;;
	v)
		VGPATH=${OPTARG}	# path to vg binary (or 'vg' if in PATH)
		;;
	?)
		echo "Invalid option: -${OPTARG}"
		exit 1
		;;
	esac
done


# Split the graph genomes with all samples set as ref into connected components (chromosomes) and put into subdirectory. Then list paths in each and rename with the path (error if >1).
mkdir ${OUTPUT}/chr_subgraphs_allref
$VGPATH chunk -C -t 4 -x ${OUTPUT}/${OUTPUT}.allrefs.xg -b ${OUTPUT}/chr_subgraphs_allref/chunk
for i in `ls ${OUTPUT}/chr_subgraphs_allref/*vg`
do
chrpaths=`$VGPATH paths -L -x ${i} | grep "${REFID}" | sed 's/.*\#.*\#//' `
echo $chrpaths | awk '{if ($2) {print "Multiple paths: "$i}}'
if [ "$chrpaths" = "" ]; then
   echo "New name not found 1: ${i}"
   exit 1;
fi
chrpaths=`echo $chrpaths | awk '{print $1}'`
if [ "$chrpaths" = "" ]; then
   echo "New name not found 2: ${i}"
   exit 1;
fi
orig=`basename $i`
echo "$orig $chrpaths" >> ${OUTPUT}/chr_subgraphs_allref/renaming.txt
mv ${i} ${OUTPUT}/chr_subgraphs_allref/${chrpaths}.vg
done

# Extract the haplotype paths (/walks) from the pangenome
IFS=',' read -r -a samplenames <<< "$NAMES"
for s in "${samplenames[@]}"
do
mkdir -p ${OUTPUT}/walks/${s}
done

mapfile -t keepchr <$KEEP
for chr in "${keepchr[@]}"
do
${VGPATH} view ${OUTPUT}/chr_subgraphs_allref/${chr}.vg | grep "^W" > ${OUTPUT}/walks/${chr}.temp.txt
for s in "${samplenames[@]}"
do
grep "${s}" ${OUTPUT}/walks/${chr}.temp.txt > ${OUTPUT}/walks/${s}/${chr}.${s}.walks.txt
done
done