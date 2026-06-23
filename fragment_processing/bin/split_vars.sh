#!/usr/bin/env bash

set -euo pipefail

# Take the counts tables annotated with variant details (concatenating multiple batches if given), sort by chromosome and first node, and split to get one line per variant (skipping fragments with no variants)
# Sorting will not be perfect but just to get an approximate order along the ref in case of batching muddling things up

# Output reorders and renames columns slightly to be consistent with previous outputs of the pipeline
# Chr	Barcode	Start	End	First_node	Strand	MAPQs	Haplotype	VarID	VarAllele	VarPos	VarType	VarCav	iPCR	SuRE_replicate1	SuRE_replicate2...



# use local tmpdir in case node runs out of space
mkdir ${PWD}/tmp
export TMPDIR=${PWD}/tmp

GAWK=gawk

OPTIND=1
OPTSTRING="o:g:"
while getopts ${OPTSTRING} opt; do
	case $opt in
	o)
		OUTNAME=${OPTARG}
		;;
	g)
		GAWK=${OPTARG}
		;;
	?)
		echo "Invalid option: -${OPTARG}"
		exit 1
		;;
	esac
done

shift $(( OPTIND - 1 ))
LISTFILES=$*

BLANKS=blanks.${OUTNAME}	# record variants where there is a blank i.e. problem getting the details in the previous step
touch ${BLANKS}

${GAWK} '
FNR == 1 {
	if (NR == 1)	{
		# print header line
		printf "Chr	Barcode	Start	End	First_node	Strand	MAPQs	Haplotype	VarID	VarAllele	VarPos	VarType	VarCav	iPCR"
		split($0,row)
		for (i=15; i<=NF; i++)	{
			printf "\t" $i
		}
		printf "\n"
	}
	next
}
{
	if ($COLVID != "NA")	{
		print $0 | "sort -k2,2 -k5,5"
	}
}
END {
    st = close("sort -k2,2 -k5,5", "to")
    if (st != 0)
        exit st
}
' ${LISTFILES} | \

${GAWK} -v blankout="${BLANKS}" '
BEGIN	{
	FS = "\t"
	OFS = "\t"
	COLVID = 9
	COLVALL = 10
	COLVPOS = 11
	COLVTYPE = 12
	COLCAV = 13
}
FNR == 1 {
	if (NR == 1)	{
		print $0
	}
	next
}
{
	if ($COLVID == "NA")	{
		next
	}
	out1 = $2 OFS $1 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS
	out2 = ""
	for (i = 14; i <= NF; i ++)	{
		out2 = out2 OFS $i
	}
	
	split($COLVID, ids, ",")
	split($COLVALL, alls, ",")
	split($COLVPOS, pos, ",")
	split($COLVTYPE, types, ",")
	split($COLCAV, cavs, ",")
	
	if ((length(ids) != length(alls)) || (length(ids) != length(pos)) || (length(ids) != length(types)) || (length(ids) != length(cavs)))	{
		print "Problem1 with lenght of variant fields: "$0 > "/dev/stderr"
		next
	}
	for (i in ids)	{
		if ((ids[i] == "") || (alls[i] == "") || (pos[i] == "") || (types[i] == "") || (cavs[i] == "")) {
			print ids[i] OFS $0 >> blankout
			continue
		}
		print out1 ids[i] OFS alls[i] OFS pos[i] OFS types[i] OFS cavs[i] out2
	}
}' | gzip > ${OUTNAME}.gz

