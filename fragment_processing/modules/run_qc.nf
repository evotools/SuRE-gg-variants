process FASTPQC	{
	label 'small_job'
	
	publishDir "$params.resDir/$params.outIpcrTrim"
	
	input:
	tuple val(reads_id), path(reads)
	
	output:
	path ("*fastqc.html")
	path "${reads_id}.fastp.json"
	
	script:
	"""
	fastqc --extract --delete --threads ${task.cpus} ${reads[0]}
	fastqc --extract --delete --threads ${task.cpus} ${reads[1]}

	fastp -i ${reads[0]} -I ${reads[1]} -A -j ${reads_id}.fastp.json --thread ${task.cpus}
	"""
}
