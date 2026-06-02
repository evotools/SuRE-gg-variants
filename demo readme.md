# SuRE-gg-variants demo
Demo for the SuRE MPRA data analysis pipelines in https://github.com/evotools/SuRE-gg-variants.git.

This demo processes SuRE data from a region of about 1 MB on bovine chromosome 3, assayed in bovine and human primary aortic endothelial cells, and identifies variants with a regulatory impact in the tested cells.

All data and config files needed for the demo run will be available at https://doi.org/10.5281/zenodo.20431190 (data archive to be released on manuscript acceptance)


## Requirements
- The demo has been tested on Linux. On Mac you may need to install the R packages for the variant analysis pipeline separately. It may run on Windows.
- ~32 GB memory required.
- To run the full demo, the cactus-pangenome command must be in the $PATH. Installing cactus can be a bit of a pain, so I have supplied the output of the cactus step so this can be skipped if inconvenient. However vg is still required, and to use with the pre-run cactus data it must be v1.61.0. A pre-compiled binary can be downloaded from https://github.com/vgteam/vg/releases/tag/v1.61.0


### Software dependencies (version tested)
- [cactus](https://github.com/ComparativeGenomicsToolkit/cactus) (2.9.3) 
- [vg](https://github.com/vgteam/vg) (1.61.0)
- Nextflow (22.10.7)
- bcftools (1.22)

## Demo Instructions
### Setup
Before running, edit the following files in the params subdirectory so they match your setup:
- [BOAEC-D1/Bi-83-A].params.yaml
	- vgPath to point to where the vg binary is installed (if in PATH just use 'vg')
	- genomePath to replace <demo_location> with the absolute path of the demo package directory
- params.analysis.demo.yaml
	- r_conf to be the absolute path of the config.cattle.demo.R file in the params subdirectory
	
Clone the github repo into the base directory
```
git clone https://github.com/evotools/SuRE-gg-variants.git
```

Most steps will be run within this subdirectory.


### Run demo
In the base demo directory, run the demo script:
```
./run_demo.sh
```

If cactus isn't straightforward, ensure you have vg v1.61.0 and run the alternative script:
```
./run_demo.no_cactus.sh
```

The expected runtime is about 20-30 minutes.

This runs the following steps:

1. Generate genome graph and associated files

Generates a genome graph of haplotype-resolved genome assemblies from the samples used to generate the SuRE libraries, plus a suitable reference genome, and produces some additional files (chromosome subgraphs, files giving individual haplotype paths in the graph) that are needed for the main pipeline.

Expected runtime ~5 minutes

2. Fragment processing pipeline

Map the reads from the SuRE libraries ('iPCR') and combine with transcript counts from cDNA sequencing files, to make a table showing fragment locations, transcript counts and associated variants.

Most time is taken setting up the conda environment (for the first run only). After that each library should be processed in a few minutes.

3. Allele-specific expression testing pipeline

For the variants in the chr3 region, test for a difference in expression betwen fragments carrying each allele.

Expected runtime 10-15 minutes including conda environment setup.


### Expected final output
Many files are generated but the main output is are three tables of variants and their allele-specific expression results, for biallelic SNPs, biallelic non-SNPs and multiallelic variants.

These are found in SuRE-gg-variants/variant_analysis/results/final_stats
