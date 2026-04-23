// Process to split processed read files into several parts to speed up gam processing
// Rename seqkit output to keep _1/_2 naming
process SPLIT_FASTQ	{
	label 'small_job'
	
	input:
	tuple val(reads_id), path(reads)
	
	output:
	path("batch*_1.fq.gz"), emit: batch_1
	path("batch*_2.fq.gz"), emit: batch_2
	val(true), emit: done
	val('trim'), emit: dir
	
	script:
	"""
	seqkit split2 -1 ${reads[0]} -2 ${reads[1]} -O ./ -p ${params.splitReads} --threads ${task.cpus}
	if (( ${params.splitReads} > 9 ))
	then
	echo "Renaming split fastq doesn't work if splitting into more than 9 parts"
	exit 1
	fi
	
	basename=`echo ${reads[0]} | sed 's/_1.fq.gz//'`
	for ((i=1;i<=${params.splitReads};i++)); do
	if [ ! -f \${basename}_1.part_00\$i.fq.gz ] || [ ! -f \${basename}_2.part_00\$i.fq.gz ]
	then
	echo "Problem splitting fastq"
	echo \${basename}_1.part_00\$i.fq.gz
	echo \${basename}_2.part_00\$i.fq.gz
	exit 1
	fi
	mv \${basename}_1.part_00\$i.fq.gz batch\$i.\${basename}_1.fq.gz
	mv \${basename}_2.part_00\$i.fq.gz batch\$i.\${basename}_2.fq.gz
	done
	
	"""
}