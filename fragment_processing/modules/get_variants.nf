// Processes to find the variants associated with each fragment using the information in the counts tables

// Find variants in the sequenced parts of the fragments
process FIND_VARS	{
	label 'single_cpu_job'
	
	publishDir "$params.resDir/$params.refName/$params.outBatchStats", pattern: "*.with_vars.${tablename}.txt"
	
	input:
	tuple val (chrid), val (tablename), path(countstable), path(variants), path(multi)
	
	output:
	tuple val (chrid), val (tablename), path("with_vars.${tablename}.txt"), path(variants), path(multi), emit: annot_countstable
	path("*.with_vars.${tablename}.txt")
	
	script:
	"""
	let gtcol=${params.gtnumber}+6
	gzip -cd ${countstable} | vars_from_countstable.pl $variants $multi \$gtcol "with_vars.${tablename}.txt"
	"""
}

process IMPUTE_VARS	{
	label 'single_cpu_job'
	
	publishDir "$params.resDir/$params.refName/$params.outBatchStats", pattern: "fragstats.*"
	
	input:
	tuple val (chrid), val (batchname), path(fragbatch), path(variants), path(multi), path(walks)
	
	output:
	tuple val (chrid), path ("imputed.${batchname}.txt"), emit: imputed_countstable
	path ("fragstats.*")
	
	script:
	"""
	let gtcol=${params.gtnumber}+6
	impute_from_countstable.pl ${walks} ${fragbatch} ${variants} ${multi} \$gtcol "imputed.${batchname}.txt"
	"""
}


// Sort fragments with variants by first node (or other sort pattern) then split into one line per variant
process SPLIT_VARS	{
	label 'single_cpu_job'
	
	publishDir "$params.resDir/$params.refName/$params.outVarCounts", pattern: "*var_counts_table.txt.gz"
	publishDir "$params.resDir/$params.refName/$params.outBatchStats", pattern: "blanks*txt"
	
	input:
	tuple val (chrid), path (fragbatches)
	
	output:
	path ("${params.libID}.${chrid}.var_counts_table.txt.gz")
	path ("blanks*txt")
	
	script:
	"""
	split_vars.sh -o ${params.libID}.${chrid}.var_counts_table.txt ${fragbatches}
	"""
	
	
}

