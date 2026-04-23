/*
Align short reads either to a linear reference genome or to a graph genome
*/

process ALIGN_LINEAR	{
	label 'big_job'
	label 'bowtie'
	
	publishDir "$params.resDir/$params.refName/$params.outAln"
		
	input:
	tuple val(reads_id), path(tr_reads)
	
	output:
	tuple val(reads_id), path("${reads_id}.namesorted.bam")
	
	script:
	"""
	bowtie2 --threads ${task.cpus} -x $params.bowtieIndex -X $params.maxInsSize -1 ${tr_reads[0]} -2 ${tr_reads[1]} | samtools view -f 2 -u | samtools sort -n -o ${reads_id}.namesorted.bam -@ ${task.cpus}
	"""
}

// Process to align paired-end reads and report alignment statistics
// Separated from the sorting step so it doesn't run out of time on the <48h nodes. This means the unsorted gam will be kept after sorting is finished and will need to be deleted if space is an issue
process GIRAFFE	{
	label 'giraffe'
	
	publishDir "$params.resDir/$params.refName/$params.outAln", pattern: 'giraffe_stats*'
	
	input:
	tuple val(reads_id), path(tr_reads)
	
	output:
	tuple val (reads_id), path("${reads_id}.unsorted.gam"), emit: alignments
	path ("giraffe_stats.${reads_id}.txt")
	
	script:
	"""
	${params.vgPath} giraffe --progress --gbz-name ${params.genomePath}.gbz --minimizer-name ${params.genomePath}.min --dist-name ${params.genomePath}.readonly.dist --fastq-in ${tr_reads[0]} --fastq-in ${tr_reads[1]} --threads ${task.cpus} > ${reads_id}.unsorted.gam
	
	${params.vgPath} stats --alignments ${reads_id}.unsorted.gam > giraffe_stats.${reads_id}.txt
	
	"""
}

// Process to sort the giraffe alignments for indexing and lookup
// Needs to run on scratch otherwise gamsort is likely to run out of space for temporary files. Could separate this step out but then need to explictly delete the unsorted file.
// Add the index to the output -- it won't actually be used in a command but needs to be in the output and the input of the next commands so that vg can find it

process GAM_SORT	{
	label 'giraffe'
	
	time { 48.hour + 24.hour * (task.attempt - 1) }
	errorStrategy 'retry'
	maxRetries 3

	
	publishDir "$params.resDir/$params.refName/$params.outAln"
	
	scratch = false
	
	input:
	tuple val (reads_id), path(unsorted_gam)
	
	output:
	tuple val (reads_id), path("${reads_id}.gam"), path("${reads_id}.gam.gai")
	
	// Make dedicated tmp directory in wd to avoid running out of space on node
	beforeScript 'export TMPDIR=$(mktemp -d --tmpdir=.)'
	//afterScript 'rm -fr $TMPDIR' uncomment when checked working as intended
	
	script:
	"""
	echo "Sorting ${unsorted_gam} for ${reads_id}"
	${params.vgPath} gamsort --index ${reads_id}.gam.gai --threads ${task.cpus} --progress ${unsorted_gam} > ${reads_id}.gam
	"""
}