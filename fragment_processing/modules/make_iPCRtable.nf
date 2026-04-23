process IPCR_BAM	{
	label 'highmem_job'
	label 'bowtie' // for using samtools
	
	
	input:
	tuple val(reads_id), path(bamfile), path(infofile)
	
	output:
	path "${reads_id}.iPCR.txt.gz"
	
	script:
	"""
	bam_to_iPCRtable.sh -b ${bamfile} -i ${infofile} -o ${reads_id}
	"""
}


// Process gam file to get iPCR table for one chromosome
// First extracts (properly paired) alignments to that chromosome from the gam file then converts that to gaf and processes with gaf_to_iPCRtable.sh.
// Input channel contains chr name, subgraph path, gam file to split, gam index (not needed for this process), info file
// The sort command as the last step before the gaf file caused problems. Currently it assumes that the read id ends in /1 or /2 and there are no /s in the previous part of the ID, uses / as the field delimiter, and sorts on field 1 (the part of the read ID common to F and R) and the first character of field 2 (1 or 2).
// Don't think this should need much processing/memory but increase on retries if failed
process IPCR_GAM	{
//	label 'highmem_job'
	cpus { 8 + (task.attempt - 1) * 4 }
	memory { task.attempt * 8.GB }
	errorStrategy 'retry'
	maxRetries 2

	publishDir "$params.resDir/$params.refName/$params.outBatchStats", pattern: "*.txt"
	
	input:
	tuple val(chr_id), path(chr_subgraph), val(reads_id), path(gamfile), path(gamindex), path(infofile)
	
	output:
	tuple val(chr_id), path("${chr_id}.${reads_id}.iPCR.txt.gz"), emit: onetable
	path ("${chr_id}.${reads_id}.skipped_pairs.txt")
	path ("${chr_id}.${reads_id}.stats.txt")
	path ("${chr_id}.${reads_id}.header")
//	path ("${chr_id}.${reads_id}.gaf") // for testing only
	
	script:
	"""
	$params.vgPath find -l ${gamfile} -A ${chr_subgraph} | \
	$params.vgPath filter --proper-pairs - | \
	$params.vgPath convert -t ${task.cpus} -G - ${chr_subgraph} | \
	sort -t/ -k1,1 -k2.1,2.1n > ${chr_id}.${reads_id}.gaf
	
	gaf_to_iPCRtable.sh -c ${chr_id} -f ${chr_id}.${reads_id}.gaf -i ${infofile} -o ${chr_id}.${reads_id}
	# rm ${chr_id}.${reads_id}.gaf
	"""
}

// Merge together iPCR tables from multiple read batches
// Output is a merged iPCR table for the chromosome, a log file of barcodes found in >1 read and what happened to them, and a list of barcodes 
// that were found in >1 read, didn't have exactly matching read alignment paths for all of them, but were kept anyway because they had the same
// nodes and fragment start/end positions
process MERGE_IPCR	{
//	label 'highmem_job'
	cpus { task.attempt * 2 }
	memory { task.attempt * 8.GB }
	errorStrategy 'retry'
	maxRetries 2
	
	publishDir "$params.resDir/$params.refName/$params.outIpcrTable", pattern: "*merged.txt.gz"
	publishDir "$params.resDir/$params.refName/$params.outBatchStats", pattern: "merge-bedpe-log.txt.gz", saveAs: { filename -> "$chr_id.$filename" }
	publishDir "$params.resDir/$params.refName/$params.outBatchStats", pattern: "multi-end-pos.txt", saveAs: { filename -> "$chr_id.$filename" }
	
	input:
	tuple val(chr_id), path ("*iPCR.txt.gz")
	
	output:
	tuple val(chr_id), path("${params.libID}.${chr_id}.merged.txt.gz"), emit: tables
	path ("merge-bedpe-log.txt.gz")
	path ("multi-end-pos.txt")
	
	script:
	"""
	merge_iPCRtables.sh -o ${params.libID}.${chr_id}.merged.txt -c ${task.cpus} *iPCR.txt.gz
	gzip merge-bedpe-log.txt
	"""
}