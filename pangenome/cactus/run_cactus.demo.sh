#!/bin/sh

# Settings for running on cluster
# #$ -N cactus
# #$ -cwd
# #$ -l h_rt=48:00:00
# #$ -l h_vmem=16G
# 
# # Need chr_list in the directory
# 
# # Initialise the environment modules
# . /etc/profile.d/modules.sh

SKIPCAC=""

MAINDIR=$1		# relative path to main github dir
SKIPCAC=$2	# if given then skip the cactus step and use the graph files in cattle-Tann-pg


# If not using the supplied graph, edit for your setup so that cactus-pangenome is in $PATH


# Needed to run
seqfile="seqfile_demo.txt"			# File giving paths to haplotype assemblies. See readme
outname="cattle-Tann-pg"				# Prefix for output names
ref="ARS_UCD2_Tann"						# Name of reference genome
allsam="ARS_UCD2_Tann Bi_83_A BOAEC_D1"	# All sample names
contigline="--refContigs 3 --otherContig chrOther"		# see readme
vgpath="vg"
	# Path to vg binary (or 'vg' if in PATH)
nonref="Bi_83_A,BOAEC_D1"	# names of the SuRE samples (not haplotypes)
chrlist="chr_list.demo.txt"		# list of contigs to include in SuRE analysis
sample_count=2				# number of samples
var_outname="var_cattle"	# prefix for variants files


# Run cactus specifying the chromosomes including X and Y and keeping other contigs
${MAINDIR}/pangenome/cactus/cactus_pangenome.demo.sh -s "$seqfile" -o "$outname" -r "$ref" -a "$allsam" -c "$contigline" -v "$vgpath" -k "$SKIPCAC" 2> cactus_out.txt
 
${MAINDIR}/pangenome/cactus/extract_walks.sh -n "$nonref" -k "$chrlist" -r "$ref" -o "$outname"

# Reformat the variants
mkdir $outname/vars
cd $outname/vars
../../${MAINDIR}/pangenome/get_var/get_variants.sh -s "$sample_count"  -n "$nonref" -v "../${outname}.vcf.gz" -o "${var_outname}" -k "../../$chrlist" -p "../../${MAINDIR}/pangenome/get_var/" -g "$vgpath"

