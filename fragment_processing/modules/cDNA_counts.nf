/*
Concatenate info files from multiple sequencing lanes and process to produce table of number of transcripts for each barcode
*/

process COUNTBARCODE	{
	label 'highmem_job'
	
	publishDir "$params.resDir/$params.outCdna", pattern: "*txt.gz"
	publishDir "$params.resDir/$params.refName/$params.outBatchStats", pattern: "*.txt"
	
	input:
	tuple val(sample_id), path (infofiles)
	
	output:
	tuple val(sample_id), path ("${sample_id}.cDNA_counts.txt.gz"), emit: counts
	path("*txt")

	script:
	"""
	count_cDNA_barcodes.sh -o ${sample_id}.cDNA_counts.txt -c ${task.cpus} $infofiles
	"""
}
