// Processes for additional processing of data after surjecting gam files to bam

// get_bedgraph.sh also includes bedGraphToBigWig
process BIGWIG	{
	label 'highmem_job'
	label 'bedgraph'
	
	publishDir "$params.resDir"
	
	input:
	tuple val(chr_id), path (counts_table)
	
	output:
	tuple val(chr_id), path("*bedGraph")
	tuple val(chr_id), path("*bw")
	
	script:
	"""
	get_bedgraph.sh -t ${counts_table} -c ${chr_id} -b ${params.bcPath} -s ${params.chrSizes}
	"""
}