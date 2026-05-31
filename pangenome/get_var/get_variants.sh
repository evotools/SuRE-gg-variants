#!/usr/bin/env bash

# get_variants.sh

VGPATH="vg"

OPTSTRING="s:n:v:o:k:g:p:"
while getopts ${OPTSTRING} opt; do
	case $opt in
	s)
		NSAMPLES=${OPTARG}	# number of samples
		;;
	n)	
		NAMES=${OPTARG}		# comma-separated list of sample names **in same order as in vcf**
		;;
	v)
		VCF=${OPTARG}		# path to cactus vcf
		;;
	o)
		OUTPUT=${OPTARG}	# output basename for cwd
		;;
	k)
		KEEP=${OPTARG}		# list of names of chromosomes to keep
		;;
	p)
		SCRIPTDIR=${OPTARG}	# path to directory with run scripts in
		;;
	g)
		VGPATH=${OPTARG}	# path to vg binary (or 'vg' if in PATH)
		;;
	?)
		echo "Invalid option: -${OPTARG}"
		exit 1
		;;
	esac
done

# Doing in two separate steps to enable checking
# Biallelic variants
# Extract fields from the vcf
let maxac=$NSAMPLES*2-1
bcftools view -g ^miss -m2 -M2 --max-ac ${maxac} ${VCF} | bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AT\t[%GT\t]\n' > vcf_extract.${OUTPUT}.txt

# Do the actual parsing
cat vcf_extract.${OUTPUT}.txt | ${SCRIPTDIR}/parse_var_path.pl biallelic.${OUTPUT} ${NAMES} y

# Multiallelic variants
bcftools view -g ^miss -Ou -a ${VCF} | bcftools view --min-alleles 3 --max-ac ${maxac}:major -Ou | bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AT\t[%GT\t]\n' > multi.vcf_extract.${OUTPUT}.txt
${SCRIPTDIR}/parse_var_path.multiallelic.pl multi.vcf_extract.${OUTPUT}.txt 142 ${NAMES} multi.${OUTPUT}

# Split up by chromosome and concatenate those from chr in the 'to keep' list (with no header in output)
mapfile -t keepchr <$KEEP
for i in biallelic multi
do
mkdir ${i}_bychr
gawk -v keep="${keepchr[*]}" -v type="${i}"	-v output="${OUTPUT}" '
BEGIN {
	split(keep, k)
	for (key in k) {
		chrkeep[k[key]]=1
	}
}
{
	if ($1 in chrkeep)	{
		print $0 >> type"_bychr/"$1"."type".txt"
		print $0 >> output"."type".chr_to_keep.txt"
	}
	else	{
		print $0 >> type".other_chr.txt"
	}
}' variants.${i}.${OUTPUT}.txt
done

