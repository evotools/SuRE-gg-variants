// Processes to combine libraries and (separately) normalise data sets

// Take vcf and extract variant IDs, unphased genotypes and allele lengths
process REFORMAT_VARS {
	label 'single_cpu_job'

	input:
	tuple val(chrname), path(vcf), path(vcf_index)

	output:
	tuple val(chrname), path("${chrname}.${params.analysis_name}.var_details.tsv")

	publishDir "$params.resdir/variant_details/", mode: 'copy'

	script:
	"""
	get_vardetails.sh -v ${vcf} -c ${chrname} \
	-o "${chrname}.${params.analysis_name}.var_details.tsv" \
	-s "${params.libHeaders}"
	"""
}

process COMBINE {
	label 'medium_job_increasing'

	input:
	tuple val(chrname), path(genotypes), path(vartables), path(mmaps)
	path(rconf)

	output:
	tuple val(chrname), path("${chrname}.${params.analysis_name}.comb.rds"), emit: comb
	path("${chrname}.${params.analysis_name}.multi.comb.rds"), emit: comb_multi
	path("${chrname}.${params.analysis_name}.counts.rds"), emit: counts
	path("${chrname}.${params.analysis_name}.vartypes.tsv"), emit: vartypes
	path("${chrname}.dup_counts.tsv.gz"), emit: mm_counts

	publishDir "$params.resdir/combined_libs/"

	script:
	"""
	set -e
	combine_libs.R ${rconf} ${chrname} ${params.analysis_name} ${genotypes} ${task.cpus}
	"""
}

process THRESHOLD {
	cpus 1
	memory 24.GB
	
	input:
	path(counts)
	path(rconf)

	output:
	path("${params.analysis_name}.thresholds.*.tsv"), emit: thresholds
	path("${params.analysis_name}.by_chr_thresholds.*.tsv")
	path("${params.analysis_name}.threshold_plots.pdf")

	publishDir "$params.resdir/filtering/", mode: 'copy'

	script:
	"""
	get_thresholds.nf.R ${rconf} ${params.analysis_name}
	"""
}

process NORMT {
	label 'single_cpu_job'

	input:
	path(lib_sizes)
	path(rconf)
	
	output:
	path("${params.analysis_name}.cDNA_normfactors.tsv")
	
	publishDir "$params.resdir/int_files/", mode: 'copy'
	
	script:
	"""
	calc_normfactors.R ${rconf} ${params.analysis_name}
	"""
}

process NORMSCORES {
	label 'normdata'

	input:
	tuple val(chrname), path(chunks)
	path(vardetails)
	path(normfactors)
	
	path(rconf)
	
	output:
	path("${chrname}.${params.analysis_name}.meanScores.tsv.gz"), emit: meanscores

	publishDir "$params.resdir/meanScores/", mode: 'copy'

	script:
	"""
	normalise_scores.R ${rconf} ${chrname} ${normfactors} ${params.analysis_name} ${task.cpus}
	"""
}

process NORMSCORES_MULTI {
	label 'normdata'

	input:
	path(chunks)
	path(normfactors)
	path(rconf)
	
	output:
	path("multi.${params.analysis_name}.meanScores.tsv.gz"), emit: meanscores_multi

	publishDir "$params.resdir/meanScores/", mode: 'copy'

	script:
	"""
	normalise_scores.multi.R ${rconf} ${normfactors} ${params.analysis_name} ${task.cpus}
	"""
}