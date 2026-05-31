#### Extract node information of variants from cactus vcf ####
## Bash script get_variants.sh to parse vcf and run scripts to classify variants.

# Use submission script to run get_variants.sh. Link appropriate chromosome list as chr_list.txt
# Edit to link to correct path for get_variants.sh
run_get_variants.cattle.sh cattle-Tann-pg.vcf.gz
run_get_variants.human.sh human-pg.vcf.gz


## Key lines in get_variants.sh
# Extract biallelic sites excluding any where all samples are non-reference (i.e. no variation within samples)
bcftools view -g ^miss -m2 -M2 --max-ac ${maxac} ${VCF} | bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AT\t[%GT\t]\n' > vcf_extract.${OUTPUT}.txt

# Run perl script to classify and get node details for biallelic variants
parse_var_path.pl ${OUTPUT}.txt ${NAMES} y

# Extract non-biallelic sites with variation within samples (N.B. May only be two alleles within the samples if e.g. 1|1 and 2|2 since ref is always counted)
bcftools view -g ^miss -Ou -a ${VCF} | bcftools view --min-alleles 3 --max-ac ${maxac}:major -Ou | bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AT\t[%GT\t]\n' > multi.vcf_extract.${OUTPUT}.txt

# Run perl script to find unique paths in each allele of the non-biallelic variants, limiting path length to 142 nodes which is the maximum possible sequenced length in iPCR
${SCRIPTDIR}/parse_var_path.multiallelic.pl multi.vcf_extract.${OUTPUT}.txt 142 ${NAMES} multi.${OUTPUT}

# Last few lines are to split up each variants file by chromosome and separate out ones on contigs you don't want to analyse.

