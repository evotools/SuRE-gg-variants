# SuRE-gg-variants
Genome graph-based analysis of regulatory variants in SuRE MPRA data

This is the repository for code used to produce the results in the cattle SuRE paper, "A high-resolution atlas of cattle regulatory variants and their cross-species activity in matched human cells", on [BioRXiv](https://doi.org/10.64898/2026.05.06.723179).

The repository has three main subdirectories:
- _pangenome_ contains code and instructions to generate the genome graph ('cactus') and extract variant lists ('get_var').
- _fragment\_processing_ contains a workflow to process one SuRE library to a table of fragment locations, transcript counts and associated variants. Its main inputs are the raw SuRE sequencing files for iPCR and cDNA, and the variant lists from the genome graph.
- _variant\_analysis_ contains a workflow to take tables from multiple SuRE libraries, combine them and test for allele-specific expression. Its subdirectory 'extra_code' contains some stand-alone scripts and code for additional processing of the results, e.g. filtering outputs and preparing plots. These include a script to produce by-chromosome lists of barcodes that map to multiple locations or appear in more than one library, which are used to exclude them by scripts in the pipeline.
- _emVar_browser_ contains an IGV session file that will load data files from the Zenodo repository to create a browsable set of variants and their SuRE results, and a readme with details of the tracks.


## Setup
### System requirements
The pipeline is designed to be run on an HPC cluster and was tested on the Edinburgh University Research Compute Cluster, using the Altair Gridengine batch system on Rocky Linux 9.

### Software dependencies (version code has been tested with)
- [cactus](https://github.com/ComparativeGenomicsToolkit/cactus) (2.9.3) 
- [vg](https://github.com/vgteam/vg) (1.61.0)
- Nextflow (22.10.7)
- bcftools (1.22)
- Anaconda

Some other tools are used but are installed in a conda environment by Nextflow when required.

### Installation
Clone the GitHub respository 
```
git clone https://github.com/evotools/SuRE-gg-variants.git
```

If dependencies above are installed then only cloning the GitHub repository is required. The first time the pipelines are run the required conda environments will be set up, taking ~10 minutes.


## To run analysis
Some adaptation to your data is required, mainly editing config and path files, and this is detailed below.

### Preparation of genome graph and variant data
Located in the pangenome subdirectory.

This step is run once for all libraries. It generates a genome graph of haplotype-resolved genome assemblies from the samples used to generate the SuRE libraries, plus a suitable reference genome, and produces some additional files (chromosome subgraphs, files giving individual haplotype paths in the graph) that are needed for the main pipeline.

Edit genome_graph/cactus/run_cactus.sh to point to the appropriate files (detailed instructions within the script and genome_graph/cactus/readme.cactus_genomes.general.txt), and to run on your compute setup.

Submit to cluster. For example to run on the UoE compute cluster
```
qsub -pe sharedmem 16 -R y run_cactus.sh
```

### Fragment processing
This is a Nextflow workflow that should be run once on each SuRE library, with a separate config file for each. Located in the fragment_processing subdirectory.

If used with an appropriate config file (e.g. conf/eddie.config) it will submit jobs to a HPC cluster, with CPUs and memory varying between jobs, mainly specified by process labels. The nextflow.config file gives reasonable settings in terms of memory and CPUs for the labels used, can be adjusted if necessary.

Parameters in conf/params.config are those that do not vary between libraries, though they may between studies (e.g. F and R adapter sequences for trimming reads), and these are loaded automatically via the nextflow.config file. Library- and/or location-specific parameters are given in a separate file supplied to nextflow via -params-file. Details of these are given in params/all_params_explained.yaml (with examples in params/Bi-83-A.yaml and params/BOAEC-D1.yaml).

Some external tools are required and are loaded via Anaconda when the conda profile is specified (see conf/conda.config).

The run is started with a single command to Nextflow, starting in the fragment_processing directory. For example to run on the UoE cluster for library Bi-83-A
```
nextflow run main.nf -c conf/eddie.config -profile conda -params-file params/Bi-83-A.yaml
```

#### Main outputs
- Trimmed reads from sequencing of SuRE plasmid libraries, and associated cutadapt 'info' files containing plasmid barcodes
- Sorted alignments of these reads with giraffe to the graph genome generated in step 1 (gam format)
- cDNA counts files -- for each barcode, the number of times it appears in the cDNA sequencing files generated from the SuRE transfections (i.e. number of transcripts)
- Counts table files, one per chromosome. Each row gives the location (with reference to the graph genome) of the gDNA fragment associated with a specific barcode, and the cDNA count for this barcode in each transfection experiment.
- Var counts tables files, one per chromosome. These are generated by finding the variants from the graph genome carried by each fragment in the counts table. One row per variant-barcode combination. This is the output required for the second pipeline.

#### Data management
Nextflow has a useful -resume flag to restart failed runs, but it requires all upstream outputs to be in place so it knows where to restart. This pipeline can generate some enormous intermediate files, particularly at the stage of aligning reads to the genome graph. Therefore it has some checkpoint options to enable starting a new run at an intermediate point in the complete pipeline, without needing to keep all intermediate outputs. Those that are required to continue the pipeline, which Nextflow would usually expect to find in a work directory, are stored in an alternative location specified by the storeDir parameter. This is required to have the same directory structure as the results directory given by resDir.

These checkpoints are, in order:
- Trim iPCR reads. If *doiPCRTrim* or *doiPCRAln* in the params file are set to "false", the cutadapt trimming step will not be rerun. Requires trimmed.\*fq.gz files and info.\*txt.gz files.
- Align iPCR reads to graph genome with giraffe. If *doiPCRAln* is set to "false" then existing \*.gam files will be used. Conversely, if *giraffeOnly* is set to "true" then the pipeline will stop after the giraffe step to allow the intermediate files to be tidied up before proceeding. In particular this generates both unsorted and sorted gam files; the former can be deleted to free up space.
- Generating counts table. If *doCountTable* is set to "false" then existing counts tables (see above) are used, and the pipeline just processes these to find the variants on each fragment.

### Variant analysis
This is a Nextflow workflow that takes the results from fragment processing for all libraries and combines them to test variants for allele-specific effects on expression. Located in the variant_analysis subdirectory.

For each variant, all fragments carrying the variant are collected and labelled by allele. Then the cDNA count for fragments carrying the variant is modelled, with allele identity included as a covariate, and a p-value for an effect of allele is reported (see publication for further details). These are collected for the whole data set and corrected for multiple testing.

#### Pre-processing
Before running through the pipeline the intermediate counts tables from all chromosomes and libraries should be run through the standalone script variant_analysis/extra_code/run_check_barcodes.both.sh to produce by-chromosome lists of barcodes that map to multiple locations or appear in more than one library, which are used to exclude them by scripts in the pipeline. Omitting this step will not cause the pipeline to fail, but you may retain data that would be better excluded.

To run this step, link all counts tables (not var counts tables) into a single directory and run the code in run_check_barcodes.both.sh. The script is set up for submission to our compute cluster and will likely need adjusted to run elsewhere. This requires a large amount of memory as it processes and sorts all the barcodes observed in the project.
```
qsub run_check_barcodes.both.sh path/to/count/table/links
```

#### Run workflow
As with the fragment processing this is set up to run on a cluster and needs a cluster config file. There are also two data-specific config files, one for Nextflow and one for R, with example files and explanations given as params/params.analysis.cattle.yaml and data/config.cattle.R, respectively. Paths to the required outputs from the fragment processing workflow are supplied to the pipeline via the Nextflow params file.

Run the pipeline in the variant_analysis subdirectory with e.g.

```
nextflow run analysis.main.nf -c conf/eddie.config -params-file params/params.analysis.cattle.yaml
```

#### Main outputs
- Tables of ASE test results, including FDR values, for all valid variants, split into SNVs, non-SNVs and multiallelic variants (results/final_stats/)
- Data files with details of fragments for significant variants (results/signif_frags/)
- Tables of various statistics for each variant, including mean SuRE scores for each allele and skew values (results/meanScores)


## Demo
Example data and full instructions for a demo data run (SuRE-gg-variants.demo.zip) can be downloaded from https://doi.org/10.5281/zenodo.20431190 (available after publication of the associated manuscript).

## Published analysis
For reference, this repository includes config files and cluster-specific scripts used to process the data in the publication above.
- pangenome/cactus/run_cactus.cattle_pg.sh
- pangenome/get_var/run_get_variants.cattle.sh
- fragment_processing/conf
- fragment_processing/data
- fragment_processing/params
- variant_analysis/conf
- variant_analyis/data
- variant_analysis/params
