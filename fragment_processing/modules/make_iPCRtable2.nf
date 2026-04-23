/*
Process to take a gam file, convert to bam file and produce iPCR table as if from the linear reference genome.
Split out into a new module so not making changes in the files for the main pipeline
*/

// Reheader is assuming contig name is in the format SAMPLE#HAP#CONTIG, renaming to CONTIG
process IPCR_GAM_BAM	{
	label 'highmem_job'
	label 'bowtie'
	
	input:
	tuple val(chr_id), path(chr_subgraph), val(reads_id), path(gamfile), path(gamindex), path(infofile)
	
	output:
	path ("*nsorted.bam"), emit: nsorted_bam
	tuple val(chr_id), path ("${chr_id}.${reads_id}.gam2bam.iPCR.txt.gz"), emit: onetable
	path ("${chr_id}.${reads_id}.gam2bam.skipped_pairs.txt")
	path ("${chr_id}.${reads_id}.gam2bam.stats")
	path ("${chr_id}.${reads_id}.gam2bam.header")
	
	script:
	"""
	$params.vgPath find -l ${gamfile} -A ${chr_subgraph} | \
	$params.vgPath filter --proper-pairs - | \
	$params.vgPath surject -x ${chr_subgraph} -b - | \
	samtools sort -@ ${task.cpus} -n -o ${chr_id}.${reads_id}.nsorted.bam
	
	gam2bam_to_iPCRtable.sh -c ${chr_id} -b ${chr_id}.${reads_id}.nsorted.bam -i ${infofile} -o ${chr_id}.${reads_id}.gam2bam
	
	if (! ${params.keepBam} )
	then
		echo "" > ${chr_id}.${reads_id}.nsorted.bam
	fi
	"""
}

// Sort bam files by position this time so they can be merged
// samtools merge will work with name-sorted files but then I will want to sort it by position afterwards anyway, so might as well do it first
process RESORT_BAM	{
	label 'medium_job'
	label 'bowtie'
	
	input:
	path(nsort_bam)
	
	output:
	path("psorted.${nsort_bam}")
	
	script:
	"""
	samtools sort -@ ${task.cpus} -o psorted.${nsort_bam} ${nsort_bam}
	"""
}

// Reheadering kept causing errors, just do manually
process MERGE_BAMS	{
	label 'medium_job'
	label 'bowtie'
	
	publishDir "$params.resDir"
	
	input:
	path(sorted_bams)
	
	output:
	path("$params.mergeBam")
	
	script:
	"""
	samtools merge --threads ${task.cpus} -o ${params.mergeBam} ${sorted_bams}
	"""
}
