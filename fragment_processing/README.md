# SuRE_SV
A Nextflow pipeline to process SuRE fragments to find their locations in the graph genome, the variants and alleles they carry, and their associated transcript counts in transfection.

The pipeline should be run once for each separate SuRE library.

Main inputs are:
- Paired-end fastq files from sequencing the iPCR library
- Single-end fastq files from sequencing the cDNAs from the transfection
- Genome graph created by cactus and its associated vcf using run_cactus scripts in cactus subdirectories.
- Walks corresponding to haplotypes from the appropriate sample extracted during the graph creation step
- Lists of variants parsed from the vcf using scripts in get_var subdirectories

Detailed input formats and additional parameters are given in params/all_params_explained.yaml.

The pipeline can be run completely from start to finish or earlier outputs can be stored elsewhere and referred to by later steps (details in all_params_explained).