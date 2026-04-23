/*
Processes to merge cDNA and iPCR tables into counts tables
*/

process MERGE_IPCR_CDNA {
//	label 'highmem_job'
	cpus { task.attempt * 2 }
	memory { task.attempt * 8.GB }
	errorStrategy 'retry'
	maxRetries 2
	
	publishDir "$params.resDir/$params.refName/$params.outCountsTables", pattern: "lib_sizes*"
	
	input:
	tuple val(chr_id), path (ipcr), val(cDNA_list)
	
	output:
	tuple val(chr_id), path ("${params.libID}.${chr_id}.counts_table_unsorted.txt.gz"), emit: table
	tuple val(chr_id), path ("lib_sizes.${params.libID}.${chr_id}.counts_table_unsorted.txt"), emit: libsizes
	
	script:
	"""
	merge_iPCR_cDNA.sh -o ${params.libID}.${chr_id}.counts_table_unsorted.txt -c ${task.cpus} -i ${ipcr} $cDNA_list
	"""
}

// Aggregate rows that align to the same place but have different barcodes
process AGG_COUNTS	{
//	label 'highmem_job'
	cpus { task.attempt * 2 }
	memory { task.attempt * 8.GB }
	errorStrategy 'retry'
	maxRetries 2
	
	publishDir "$params.resDir/$params.refName/$params.outCountsTables"
	
	input:
	tuple val(chr_id), path(countstable)
	val(samplecount)
	
	output:
	tuple val(chr_id), path ("${params.libID}.${chr_id}.counts_table.txt.gz"), emit: counts_table
	path "${params.libID}.${chr_id}.counts_table.txt_pos_multi_BC.txt.gz"
	path "${params.libID}.${chr_id}.counts_table.txt_merged*.txt.gz"
	
	script:
	"""
	CMD="aggregate_counts.sh -u $countstable -o ${params.libID}.${chr_id}.counts_table.txt -n $samplecount"
	if [ ! -z $params.countTableSort ]
	then
	  CMD="\${CMD} -s $params.countTableSort"
	fi
	eval \$CMD
	"""
}

// Find how many unique fragments and totals for iPCR and each library
process COUNT_OUTPUT	{
	label 'medium_job'
	
	publishDir "$params.resDir/$params.refName/$params.outCountsTables"
	
	input:
	tuple val(chr_id), path ("${params.libID}.${chr_id}.counts_table.txt.gz")
	
	output:
	path "${params.libID}.${chr_id}.lib_counts.txt"
		
	script:
	"""
	get_lib_sizes.R ${params.libID}.${chr_id}.counts_table.txt.gz ${params.libID}.${chr_id}.lib_counts.txt
	"""
}
