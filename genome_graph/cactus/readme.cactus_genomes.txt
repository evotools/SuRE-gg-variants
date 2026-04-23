#### Get data for graph genomes and align using cactus ####

### Human
## Download reference genome and haplotype-resolved assembly for NA12878
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
wget ftp://ftp.dfci.harvard.edu/pub/hli/whdenovo/asm/NA12878-denovo*

## Make seqfile.hum.txt according to 
# https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#interface
# Contents of file looks like:
CHM13	/exports/eddie/scratch/lplender/data/human_haps/hs1.fa.gz
NA12878.1	/exports/eddie/scratch/lplender/SuRE/data/human_haps/NA12878-denovo-H1.fa.gz
NA12878.2	/exports/eddie/scratch/lplender/SuRE/data/human_haps/NA12878-denovo-H2.fa.gz

## Submit job to run minigraph-cactus on eddie
qsub -pe sharedmem 16 -R y <path-to-SuRE_SV>/pangenome/cactus/run_cactus.human_pg.sh


### Cattle
## Download reference genome ARS-UCD2.0 with chromosomes named with numbers from BovLRC
wget https://bovlrc.s3.amazonaws.com/ARS_UCD_v2.0.fa.gz

# Download T. annulata genome and concatenate with renamed cattle genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/003/225/GCA_000003225.1_ASM322v1/GCA_000003225.1_ASM322v1_genomic.fna.gz
mv ARS_UCD_v2.0.fa.gz ARS_UCD_v2.0_Tann.fa.gz
cat GCA_000003225.1_ASM322v1_genomic.fna.gz >> ARS_UCD_v2.0_Tann.fa.gz

## Convert gfa output from hifiasm to fasta
for i in BOAEC_D1 Bi_TaA
do
for j in 1 2
do
zcat ${i}_trim.asm.bp.hap${j}.p_ctg.gfa.gz | awk '/^S/{print ">"$2;print $3}' | gzip > ${i}_trim.asm.bp.hap${j}.p_ctg.fa.gz
done
done

## Make seqfile.cattle.txt
# Sample names cannot contain . so just use ARS_UCD2_Tann as the reference name
# Bi_TaA is Bos indicus (Sahiwal), BOAEC is Bos taurus (Holstein-Friesan)
# Sample names (in seqfile are followed by paths):
ARS_UCD2_Tann
Bi_83_A.1
Bi_83_A.2
BOAEC_D1.1
BOAEC_D1.2

## Submit job to run minigraph-cactus on eddie
qsub -pe sharedmem 16 -R y <path-to-SuRE_SV>/pangenome/cactus/run_cactus.cattle_pg.sh


### Split into chromosome subgraphs
# To make processing of the alignments more manageable, need to split into chromosomes
# The way I do this is to split the graph genome into its connected components 
# (which correspond to chromosomes in the reference) and then use these in the pipeline to
# extract the reads aligned to each chromosome using vg find.
# The following command (included in cactus_pangenome.sh) splits the graph into subgraphs
# which are then renamed according the chromosome they represent (the reference path they contain).
vg chunk -C -t <threads> -b <basename> -x human-pg.gbz


### Extract walks (for later variant/fragment processing)
# Again, split into connected components but this time from the version of the pangenome with all samples set to ref
# Then vg view each of the subgraphs, parse the walks and then the specific sample walks out of the GFA produced
# Done by extract_walks.sh
