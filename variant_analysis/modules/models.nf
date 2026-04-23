// Prepare model data and run models

process PREP_DATA {
	label 'big_job_increasing'

	input:
	tuple val(chrname), path(combdata)
	path(thresholds)
	path(rconf)
	
	output:
	tuple val (chrname), path("chunk.${chrname}_*.tsv.gz"), emit: chunks
	path("${chrname}.${params.analysis_name}.filter_details.tsv")
	
	publishDir "$params.resdir/filtering/", mode: 'copy', pattern: "*filter_details.tsv"
	publishDir "$params.resdir/chunks/", pattern: "chunk.*.tsv.gz"
	
	script:
	"""
	prep_data.nf.R ${rconf} ${chrname} ${combdata} ${params.analysis_name} ${task.cpus}
	"""
}

// Run models for one species (or all), on one chunk
// Outputs model data plus a list of variants for which the model fitting failed
process RUN_MODELS {
	label 'models'
	
	input:
	tuple val(chrname), path(chunkfile), val(species)
	path(rconf)
	
	output:
	path("model_failed*tsv")
	//tuple val(species), val(chrname), path("modelResults.*tsv.gz"), emit: modelRes
	path("modelResults.*tsv.gz"), emit: modelRes
	
	publishDir "$params.resdir/model_results"
	
	script:
	"""
	run_models.pll.nf.R ${rconf} ${species} ${chrname} ${chunkfile} ${params.analysis_name}
	"""
}

process PREP_DATA_MULTI	{
	label 'big_job_increasing'

	input:
	path(multi_comb)
	path(gts_locs)
	path(thresholds)
	path(rconf)

	output:
	path("chunk_multi.*.tsv.gz"), emit: multi_chunks
	path("chunk.*.tsv.gz"), emit: extra_chunks
	path("${params.analysis_name}.multi_now_bi.recoding.tsv"), emit: recoding
	path("${params.analysis_name}.filter_details_multi.tsv")
	path("${params.analysis_name}.threshold_plots.multi.pdf")
	path("${params.analysis_name}.multi_now_bi.recoded_gt.tsv"), emit: recodegt

	publishDir "$params.resdir/filtering/", pattern: '*.{pdf,tsv}'
	publishDir "$params.resdir/chunks_multi/", pattern: "chunk_multi.*.tsv.gz"
	publishDir "$params.resdir/chunks/", pattern: "chunk.*.tsv.gz"
	
	script:
	"""
	prep_data.multi.nf.R ${rconf} ${params.analysis_name} ${task.cpus}
	"""
}

process RUN_MODELS_MULTI {
	label 'models_multi'

	input:
	tuple path(chunkfile), val(species)
	path(rconf)
	
	output:
	path("modelResults.*tsv.gz"), emit: modelRes
	path("model_failed*tsv")

	publishDir "$params.resdir/model_results"

	script:
	"""
	run_models_multi.pll.nf.R ${rconf} ${species} ${chunkfile} ${params.analysis_name}
	"""
}