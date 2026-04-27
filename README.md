# SuRE-gg-variants
Genome graph-based analysis of regulatory variants in SuRE MPRA data

This is the repository for code used to produce the results in the cattle SuRE paper, "A high-resolution atlas of cattle regulatory variants and their cross-species activity in matched human cells"

Directory 'genome_graph' contains code and instructions to generate the genome graph ('cactus') and extract variant lists ('get_var').

Directory 'fragment_processing' contains the 'SuRE_SV' pipeline to process one SuRE library to a table of fragment locations, transcript counts and associated variants. Its main inputs are the raw SuRE sequencing files for iPCR and cDNA, and the variant lists from the genome graph.

Directory 'variant_analysis' contains the 'SuRE-gg-analysis' pipeline to take tables from multiple SuRE libraries, combine them and test for allele-specific expression. Its subdirectory 'extra_code' contains some stand-alone scripts and code for additional processing of the results, e.g. filtering outputs and preparing plots.
Before running through the pipeline the intermediate counts tables from all chromosomes and libraries should be run through the standalone script variant_analysis/extra_code/run_check_barcodes.both.sh to produce by-chromosome lists of barcodes that map to multiple locations or appear in more than one library, which are used to exclude them by scripts in the pipeline.