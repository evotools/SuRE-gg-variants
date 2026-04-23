/*
Uses cutadapt on iPCR reads to
	- Trim off F and R adapters at start and end
	- Add /1 or /2 suffix to read name
	- Remove reads which are longer than they should be, i.e. expected sequences are not all present (read length minus adapter length, and for F reads minus barcode length)
	- Check for readthrough of the fragment to the F and R adapter at the other end and trim this off
	- Remove reads of less than $params.minsize bp after trimming
	- Produce info file to extract barcodes
	
	
For cDNA reads, a single round of cutadapt to
	- Trim F adapter (different from iPCR) off start of single-end read
	- Produce info file to count barcodes
	- Remove trimmed reads file as not needed?
*/


process CUTADAPT {
	label 'medium_job'
	label 'cutadapt'
	
	publishDir "$params.resDir/$params.outIpcrTrim"
	
	input:
	tuple val(reads_id), path(reads)
	
	output:
	tuple val(reads_id), path ("trimmed.${reads_id}_*.fq.gz"), emit: tr_reads
	tuple val(reads_id), path ("info.${reads_id}.txt.gz"), emit: info
	path ("stats*"), emit: reports
	
	script:
	"""
	cutadapt --cores 0 -g $params.adapterF -G $params.adapterR\\;max_error_rate=0.13 -O 4 --maximum-length 89:142 -m $params.minsize --rename='{id}/{rn}' --too-long-output too_long.${reads_id}_1.fq.gz --too-long-paired-output too_long.${reads_id}_2.fq.gz --info-file=info.${reads_id}.txt.gz -o ${reads_id}_1.flt1.fq.gz -p ${reads_id}_2.flt1.fq.gz ${reads[0]} ${reads[1]} > stats.${reads_id}.main.txt
	
 	fadaptRC=\$( echo $params.adapterF | grep '^[atcgATCG]' | rev | tr atcgATCG tagcTAGC )
 	fadaptRCint=\${fadaptRC}NNNNNNNNNNNNNNNNNNNNX
 	radaptRC=\$( echo $params.adapterR | grep '^[atcgATCG]' | rev | tr atcgATCG tagcTAGC )
 	radaptRCint=\${radaptRC}"X"
	
	cutadapt --cores 0 -a \${radaptRCint}\\;max_error_rate=0.13 -A \${fadaptRCint} -O 4 -m $params.minsize -o trimmed.${reads_id}_1.fq.gz -p trimmed.${reads_id}_2.fq.gz ${reads_id}_1.flt1.fq.gz ${reads_id}_2.flt1.fq.gz > stats.${reads_id}.3prime.txt
	"""
}

process CUTADAPT_CDNA	{
	label 'medium_job'
	label 'cutadapt'
	
	publishDir "$params.resDir/$params.refName/$params.outBatchStats", pattern: "stats_cDNAcutadapt.*"

	input:
	tuple val(sample_id), path(readfile), val(readsID)
	
	output:
	tuple val (sample_id), path ("info.${sample_id}.${readsID}.txt.gz")
	
	script:
	"""
	cutadapt --cores 0 -g $params.cDNAadapterF --info-file=info.${sample_id}.${readsID}.txt.gz -o trimmed.${readsID}.fq.gz ${readfile} > stats_cDNAcutadapt.${readsID}.txt
	"""
}
