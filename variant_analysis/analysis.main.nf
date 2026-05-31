/*
Pipeline to take the var_counts tables produced from SuRE data, prepare the fragment data
for analysis, and run models to find variants with a significant difference between alleles.
Options to start from the beginning (combining libraries), preparing the data (so can use
different filtering options), and running the models.
*/

// Needed to run on eddie
nextflow.enable.dsl=2

// Import processes and sub-workflows
include { NORMT; COMBINE; THRESHOLD; REFORMAT_VARS; NORMSCORES; NORMSCORES_MULTI } from './modules/combine_and_norm.nf'
include { PREP_DATA; RUN_MODELS; PREP_DATA_MULTI; RUN_MODELS_MULTI } from './modules/models.nf'
include { GET_ALLELES; NOWBIALL_DETAILS; GET_PLOTDATA; GET_PLOTDATA_MULTI; SIGNIF_FRAGS; SIGNIF_FRAGS_MULTI; VAR_PAIRS } from './modules/results.nf'


workflow	{
	// Get the variant details from the vcf to save parsing multiple times
	// Makes one file per chromosome on the list given by chrList
	// Assume there's a vcf index with the same path plus .tbi
	def list_of_chr = params.chrList.split(" ") as List
	Channel
  		.fromPath(params.vcfPath)
 		.flatMap { vcf ->
    		def tbi = "${vcf}.tbi"
    		list_of_chr.collect { chr -> tuple(chr, vcf, tbi) }
  		}
  	.set { chr_vcf_ch }
	REFORMAT_VARS(chr_vcf_ch)

	// Get the var counts table files
	vartable_ch = Channel.fromPath(params.vartable_paths).flatMap()
	vartable_ch
		.map { it ->
			def base = it.getBaseName().replaceFirst(/\.var_counts_table\.txt/, '')
			def chrname = base.tokenize('.')[1]  // get the second part (index 1)
			[chrname, it]
		}
		.groupTuple()
		.set { vartable_ch }
	
	// Add variant details
	REFORMAT_VARS.out
		.join(vartable_ch)
		.set { vartable_gt_ch }

	// Get the multi-mapping barcodes files and join to the vartable files
	mmap_ch = Channel.fromPath(params.mmap_path).flatMap()
	mmap_ch
		.map { it ->
			def base = it.getBaseName().replaceFirst(/\.dup_bc\.tsv/, '')
			def chrname = base.tokenize('.')[1]
			[chrname, it]
		}
		.set { mmap_ch }
	
//	mmap_ch.view()
	vartable_gt_ch = vartable_gt_ch.join(mmap_ch)
	COMBINE(vartable_gt_ch, params.r_conf)

	// Get thresholds for filtering
	counts_ch = COMBINE.out.counts.collect()
	THRESHOLD(counts_ch, params.r_conf)

	// Prepare model data - filter and split into chunks
	//COMBINE.out.norm.view()
	PREP_DATA(COMBINE.out.comb, THRESHOLD.out.thresholds, params.r_conf)

	// Separate the chunks output into separate elements then make one for each species/all
	def list_of_spec = params.speciesList.split(" ") as List
	PREP_DATA.out.chunks
		.flatMap { chr, files -> files.collect { f -> tuple(chr, f) } }
		.flatMap { chr, file -> list_of_spec.collect { spec -> tuple(chr, file, spec) } }
		.set { formodels_init_ch }
	
	// Prepare the data for non-biallelic variants
	multi_ch = COMBINE.out.comb_multi.collect()
	REFORMAT_VARS.out
		.map { chr, file -> file }
		.collect()
		.set { allgt_ch }
	PREP_DATA_MULTI(multi_ch, allgt_ch, THRESHOLD.out.thresholds, params.r_conf)

	// Add the now-biallelic chunks to the channel for biallelic models
	PREP_DATA_MULTI.out.extra_chunks
		.flatMap { it }
		.map { Path c ->
			def m = (c.name =~ /chunk\.(.*?)_\d+\.tsv.gz/)
			def chr = m[0][1]
			tuple( chr, c )
		}
		.set { extra_chunks_ch }
	extra_chunks_ch
		.flatMap { chr, file -> list_of_spec.collect { spec -> tuple(chr, file, spec) } }
		.concat( formodels_init_ch )
		.set { formodels_ch }
	//formodels_ch.view()

	// Run model tests on the chunks. Rerun failures with more resources
	RUN_MODELS(formodels_ch, params.r_conf)
	RUN_MODELS.out.modelRes
		.collect()
		.set { modelres_ch }
	
	// Set up the model channel for the non-biallelic variants and run models
	PREP_DATA_MULTI.out.multi_chunks
		.flatMap { it }
		.flatMap { file -> list_of_spec.collect { spec -> tuple(file, spec) } }
		.set { formodels_multi_ch }
	RUN_MODELS_MULTI(formodels_multi_ch, params.r_conf)
	RUN_MODELS_MULTI.out.modelRes
		.collect()
		.set { modelres_multi_ch }

	// Get the library sizes and calculate a table for normalising
	libsizes_ch = Channel.fromPath(params.libsize_paths).collect()
	NORMT(libsizes_ch, params.r_conf)

	// Find the mean normalised scores for each variant, from the chunks (filtered fragments)
	// for each chromosome
	PREP_DATA.out.chunks
		.flatMap { chromosomeID, filePathList -> filePathList.collect { filePath -> [chromosomeID, filePath] } }
		.concat(extra_chunks_ch)
		.groupTuple()
		.set { all_chunks_ch }
	// Lazily give it all vardetails and R script will select the correct one
	REFORMAT_VARS.out
		.map {it -> it[1]}
		.collect()
		.set { vardetails_ch }
	NORMSCORES(all_chunks_ch, vardetails_ch, NORMT.out, params.r_conf)

	// Get single file of normalised scores for non-biallelic variants
	NORMSCORES_MULTI(PREP_DATA_MULTI.out.multi_chunks, NORMT.out, params.r_conf)
	
	// Extract some extra details from the vcf
	GET_ALLELES(Channel.fromPath(params.vcfPath))

	// Get details of the alleles that were actually compared in the variants that were annotated as multiallelic
	// but contained only two alleles
	NOWBIALL_DETAILS(PREP_DATA_MULTI.out.recoding, GET_ALLELES.out.alleles, params.r_conf)

	// Process the model results, correct for multiple tests,
	// output with mean scores for plotting
	GET_PLOTDATA(modelres_ch, allgt_ch, NORMSCORES.out.meanscores.collect(), COMBINE.out.vartypes.collect(), 
	GET_ALLELES.out.alleles, NOWBIALL_DETAILS.out.nowbi_dat, PREP_DATA_MULTI.out.recodegt, params.r_conf)
	GET_PLOTDATA_MULTI(modelres_multi_ch, allgt_ch, GET_ALLELES.out.alleles, NORMSCORES_MULTI.out.meanscores_multi, params.r_conf)

	// For plotting, extract the fragments for each significant variant, normalise scores and get barcodes
	vartable_ch
		.concat(all_chunks_ch)
		.groupTuple()
		.flatMap { chrid, listlist -> [tuple(chrid, listlist.flatten())] }
		.set { vartable_chunks_ch }
	SIGNIF_FRAGS(vartable_chunks_ch, GET_PLOTDATA.out.snpsfile, GET_PLOTDATA.out.otherfile, NORMT.out, params.r_conf)

	SIGNIF_FRAGS_MULTI(vartable_ch, GET_PLOTDATA_MULTI.out.multifile, PREP_DATA_MULTI.out.multi_chunks, NORMT.out, params.r_conf)

	// Run the script to look at pairs of variants on the same fragment
	VAR_PAIRS(SIGNIF_FRAGS.out.signif_chunks, GET_PLOTDATA.out.snpsfile, GET_PLOTDATA.out.otherfile, params.r_conf)

}