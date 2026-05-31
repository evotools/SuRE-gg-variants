#### Get data for graph genomes and align using cactus ####

## Make seqfile.txt according to 
https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#interface
Contents of file looks like:

CHM13	/exports/eddie/scratch/lplender/data/human_haps/hs1.fa.gz
NA12878.1	/exports/eddie/scratch/lplender/SuRE/data/human_haps/NA12878-denovo-H1.fa.gz
NA12878.2	/exports/eddie/scratch/lplender/SuRE/data/human_haps/NA12878-denovo-H2.fa.gz

## Update run scripts
Copy run_cactus.sh and adapt to run with your setup and filenames

## Example command -- to submit job to run minigraph-cactus and extract walks for each haplotype, on the UoE compute cluster, with an appropriate number of cores
> qsub -pe sharedmem 16 -R y run_cactus.sh


### Tasks performed by the scripts

## Run cactus-pangenome, generating files needed for giraffe aligner
For details see https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md

## Make additional files related to the graph genome
Read-only version of distance file needed by giraffe (otherwise fails when running multiple batches at once)
Vg version of graph
Index of graph
Version of graph with all samples set to be reference

## Split into chromosome subgraphs
To make processing of the alignments more manageable, need to split into chromosomes
The way I do this is to split the graph genome into its connected components 
(which correspond to chromosomes in the reference) and then use these in the pipeline to
extract the reads aligned to each chromosome using vg find.
The following command (included in cactus_pangenome.sh) splits the graph into subgraphs
which are then renamed according the chromosome they represent (i.e. the reference contig they contain).
> vg chunk -C -t <threads> -b <basename> -x human-pg.gbz

## Extract walks (for later variant/fragment processing). Done by extract_walks.sh
Again, split into connected components but this time from the version of the pangenome with all samples set to ref. Then vg view each of the subgraphs, parse the walks and then the specific sample walks out of the GFA produced.

