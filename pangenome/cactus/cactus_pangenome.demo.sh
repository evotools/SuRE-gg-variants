#!/usr/bin/env bash

# cactus_pangenome.sh

# Some edits specifically for the demo

mkdir ${PWD}/tmp
export TMPDIR=${PWD}/tmp

REFCONTIGS=""
SKIPCACTUS=""
OPTIND=1
OPTSTRING="s:o:r:a:c:v:k:"
while getopts ${OPTSTRING} opt; do
	case $opt in
	s)
		SEQLIST=${OPTARG}
		;;
	o)
		OUTPUT=${OPTARG}  # for both output directory and name of outputs within it, in wd
		;;
	r)
		REF=${OPTARG} # single reference for main output
		;;
	a)
		ALLREFS=${OPTARG}
		;;
	c)
		REFCONTIGS=${OPTARG}
		;;
	v)
		VGPATH=${OPTARG}
		;;
	k)
		SKIPCACTUS=${OPTARG}
		;;
	?)
		echo "Invalid option: -${OPTARG}"
		exit 1
		;;
	esac
done


# Run the main graph genome construction, change memory if required
CMD="cactus-pangenome ./js ./${SEQLIST} --consMemory 8G --indexMemory 8G --outDir ${OUTPUT} --outName ${OUTPUT} --reference ${REF} --giraffe clip filter --vcf --gfa --gbz --viz ${REFCONTIGS}"
echo ${CMD}
if [ -z "$SKIPCACTUS" ]
then
	eval ${CMD}
else
	echo "Skipping cactus-pangenome call"
fi

# Make a read-only version of the distance file needed by giraffe
cp ${OUTPUT}/${OUTPUT}.dist ${OUTPUT}/${OUTPUT}.readonly.dist
chmod -w  ${OUTPUT}/${OUTPUT}.readonly.dist

# Make the graph genome available in different formats
$VGPATH convert -t 8 ${OUTPUT}/${OUTPUT}.gbz > ${OUTPUT}/${OUTPUT}.vg
$VGPATH index -x ${OUTPUT}/${OUTPUT}.xg -t 8 ${OUTPUT}/${OUTPUT}.gbz

# Make index with all as references
$VGPATH gbwt -Z --set-tag "reference_samples=${ALLREFS}" --gbz-format -g ${OUTPUT}/${OUTPUT}.allrefs.gbz ${OUTPUT}/${OUTPUT}.gbz
$VGPATH index -x ${OUTPUT}/${OUTPUT}.allrefs.xg -t 8 ${OUTPUT}/${OUTPUT}.allrefs.gbz

# Split into connected components (chromosomes) and put into subdirectory. Then list paths in each and rename with the path (error if >1).
mkdir ${OUTPUT}/chr_subgraphs
$VGPATH chunk -C -t 4 -x ${OUTPUT}/${OUTPUT}.gbz -b ${OUTPUT}/chr_subgraphs/chunk
for i in `ls ${OUTPUT}/chr_subgraphs/*vg`
do
chrpaths=`$VGPATH paths -L -x ${i} | grep ${REF} | sed 's/.*\#.*\#//'`
echo $chrpaths | awk -v filename=$i '{if ($2) {print "Multiple paths: " filename}}'
chrpaths=`echo $chrpaths | awk '{print $1}'`
orig=`basename $i`
echo "$orig $chrpaths" >> ${OUTPUT}/chr_subgraphs/renaming.txt
mv ${i} ${OUTPUT}/chr_subgraphs/${chrpaths}.vg
done
