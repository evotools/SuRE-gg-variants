# SuRE-gg-analysis
A Nextflow pipeline to test variants for allele-specific expression in SuRE.

This pipeline can be used to combine multiple libraries processed with the SuRE_SV pipeline, as long as they all used the same genome graph (i.e. tested the same set of fragments).

Before running through the pipeline the intermediate counts tables from all chromosomes and libraries should be run through the standalone script run_check_barcodes.both.sh to produce by-chromosome lists of barcodes that map to multiple locations or appear in more than one library, which are used to exclude them by scripts in the pipeline.
These files list all valid fragments and their locations, regardless of whether they carried variants or not, and are found in the counts_tables subdirectory of the SuRE_SV results directory.

The pipeline works on each chromosome to combine the fragment details from each library, collect fragments for each variant, do fragment and variant filtering, and test variants for allele-specific expression.

Its main input is the SuRE variant counts tables, found in the var_counts_tables subdirectory of the SuRE_SV results directory. These files contain one line for each variant on each fragment, and give the cDNA counts associated with the fragment from each transfection. Files giving the sizes of each library are also used; these are found in the SuRE_SV lib_sizes subdirectory.

Other inputs are described in data/config.cattle.eddie.R and params/params.analysis.cattle.yaml.