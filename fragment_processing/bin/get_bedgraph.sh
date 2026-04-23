#!/usr/bin/env bash

# get_bedgraph.sh

OPTIND=1
OPTSTRING="t:c:b:s:"
while getopts ${OPTSTRING} opt; do
	case $opt in
	t)
		TABLE=${OPTARG}
		;;
	c)
		CHR=${OPTARG}
		;;
	b)
		BED2GRAPH=${OPTARG}
		;;
	s)
		CHRSIZES=${OPTARG}
		;;	
	?)
		echo "Invalid option: -${OPTARG}"
		exit 1
		;;
	esac
done

header=$(gzip -cd ${TABLE} | head -n 1)
read -r -a header_arr <<< "$header"

for i in `seq 8 "${#header_arr[@]}"`
do
countcol=$((${i}-1))
CMD="gzip -cd ${TABLE} | grep -v \"chr\" | cut -f 2,3,4,${i} | ${BED2GRAPH} > ${CHR}.${header_arr[$countcol]}.bedGraph"
echo ${CMD2}
eval ${CMD}
CMD2="bedGraphToBigWig ${CHR}.${header_arr[$countcol]}.bedGraph ${CHRSIZES} ${CHR}.${header_arr[$countcol]}.bw"
echo ${CMD2}
eval ${CMD2}
done
