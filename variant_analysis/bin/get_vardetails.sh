#!/usr/bin/env bash

# get_vardetails.sh

# Takes vcf file and extracts ref and all alt lengths, and converts phased genotypes to
# unphased (including ordering numerically so 3|1 becomes 1/3)
# Will return information for all variants with no missing data
# If sample headers are supplied then will use these (and assume they're the same order as in vcf),
# otherwise will use sample names from the vcf

OPTIND=1
OPTSTRING="v:c:o:s:"
while getopts ${OPTSTRING} opt; do
	case $opt in
	v)
		VCF=${OPTARG}
		;;
	c)
		CHRID=${OPTARG}
		;;
	o)
		OUTNAME=${OPTARG}
		;;
	s)
		SNAMES=${OPTARG}
		;;
	?)
		echo "Invalid option: -${OPTARG}"
		exit 1
		;;
	esac
done

if [ -z "${SNAMES}" ]; then
SNAMES="$(bcftools query -l "${VCF}" | paste -sd '\t' -)"
fi

echo "VarID	Chr	RefLength	AltLengths	${SNAMES}" > ${OUTNAME}

bcftools view -g ^miss ${VCF} ${CHRID}| bcftools query -f "%ID\t%CHROM\t%REF\t%ALT[\t%GT]\n" \
	| gawk '
	{
	n = split($4,alts,",")
	altlens = length(alts[1])
	for (i = 2; i<=n; i++) {
		altlens = altlens "," length(alts[i])
	}
	printf "%s\t%d\t%d\t%s", $1, $2, length($3), altlens
	for (i = 5; i <= NF; i++) {
		split($i, alleles ,"|")
		asort(alleles)
		printf "\t%s", alleles[1]"/"alleles[2]
	}
	printf "\n"
	}' >> ${OUTNAME}
