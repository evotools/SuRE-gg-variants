#!/bin/sh
#$ -N cactus
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=16G

# Run with
# qsub -pe sharedmem 16 -R y run_cactus.sh

# Initialise the environment modules
. /etc/profile.d/modules.sh
module load roslin/bcftools/1.20

MAINDIR=$1	# absolute path to genome_graph subdirectory

# Load virtualenv for cactus
# Edit for your setup so that cactus-pangenome is in $PATH
source /exports/cmvm/eddie/eb/groups/prendergast_roslin/Lindsey/programs/cactus/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate

# Needed to run
seqfile="seqfile.cattle.txt"			# File giving paths to haplotype assemblies. See readme
outname="cattle-Tann-pg"				# Prefix for output names
ref="ARS_UCD2_Tann"						# Name of reference genome
allsam="ARS_UCD2_Tann Bi_83_A BOAEC_D1"	# All sample names
contigline="--refContigs 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 X Y CR940347.1 CR940349.1 CR940350.1 CR940351.1 CR940348.1 CR940352.2 CR940353.1 CR940346.1 --otherContig chrOther"		# see readme
vgpath="/exports/cmvm/eddie/eb/groups/prendergast_roslin/Lindsey/programs/vg/vg"
	# Path to vg binary
nonref="Bi_83_A,BOAEC_D1"	# names of the SuRE samples (not haplotypes)
chrlist="chr_list.txt"		# list of contigs to include in SuRE analysis
sample_count=2				# number of samples
var_outname="var_cattle"	# prefix for variants files


# Run cactus specifying the chromosomes including X and Y and keeping other contigs
${MAINDIR}/cactus/cactus_pangenome.any.sh -s "$seqfile" -o "$outname" -r "$ref" -a "$allsam" -c "$contigline" -v "$vgpath"
 
${MAINDIR}/cactus/extract_walks.sh -n "$nonref" -k "$chrlist" -r "$ref" -o "$outname" -v "$vgpath"

mkdir $outname/vars
cd $outname/vars
${MAINDIR}/get_var/get_variants.sh -s "$sample_count"  -n "$nonref" -v "../${outname}.vcf.gz" -o "${var_outname}" -k "../../$chrlist" -p "${MAINDIR}/get_var/" -g "$vgpath"
