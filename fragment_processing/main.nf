/*
Nextflow workflow to analyse SuRE data to look at the effects of structural variants
*/

// Needed to run on eddie
nextflow.enable.dsl=2

// Import processes and sub-workflows
include { CUTADAPT } from './modules/cutadapt.nf'
include { FASTPQC as FASTQC1 } from './modules/run_qc.nf'
include { FASTPQC as FASTQC2 } from './modules/run_qc.nf'
include { SPLIT_FASTQ } from './modules/split_fastq.nf'
include { GIRAFFE; GAM_SORT } from './modules/align_sr.nf'
include { IPCR_GAM; MERGE_IPCR } from './modules/make_iPCRtable.nf'
include { COUNTS } from './workflows/counts.nf'
include { FIND_VARS; IMPUTE_VARS; SPLIT_VARS } from './modules/get_variants.nf'
include { COUNT_OUTPUT } from './modules/counts_tables.nf'

workflow {
	// PART 1 - Align iPCR reads to graph genome
	// Get trimmed iPCR reads, either from trimming raw reads or finding reads that have already been trimmed
	if ( params.doiPCRAln )	{
		if ( params.doiPCRTrim ) {
			println "Trimming reads"
			params.readpattern = "${params.iPCRDir}/*{1,2}.fq.gz"
			ch_read_pairs = Channel.fromFilePairs(params.readpattern, checkIfExists: true)
			FASTQC1(ch_read_pairs)
		
			// If needed, split the trimmed fastq into batches to speed up processing
			// parameter splitReads gives number of batches to split into
			// New identifier used to name outputs is readID string plus batch number
			if (params.splitReads)	{
				println "Splitting reads before trim"
				SPLIT_FASTQ (ch_read_pairs)
				SPLIT_FASTQ.out.batch_1
					.flatten()
					.map { it -> [it.getBaseName().replaceFirst(/_1.fq/,""), it] }
					.set { ch_batches_1 }
				SPLIT_FASTQ.out.batch_2
					.flatten()
					.map { it -> [it.getBaseName().replaceFirst(/_2.fq/,""), it] }
					.join(ch_batches_1, failOnMismatch: true)
					.map { it -> [it[0], [it[2], it[1]]] }
					.set { ch_read_pairs }
			}
		
			CUTADAPT(ch_read_pairs)
			ch_iPCRreads = CUTADAPT.out.tr_reads
			ch_iPCRinfo = CUTADAPT.out.info
			FASTQC2(ch_iPCRreads)
		}
		else	{
		// check the storage directory for fastq files using fromFilePairs, then change the key so it matches the original read ID
			println "Finding pre-trimmed reads"
			params.trimpattern = "${params.storeDir}/${params.outIpcrTrim}/*{1,2}.fq.gz"
			
			Channel
				.fromFilePairs(params.trimpattern, checkIfExists: true)
				.map { it -> [it[0].replaceFirst(/trimmed./, ""), it[1]] }
				.set { ch_iPCRreads }
		
		// Find the info files
		params.infopattern = "${params.storeDir}/${params.outIpcrTrim}/info*txt.gz"
		Channel
			.fromPath(params.infopattern, checkIfExists: true)
			.map { it -> [it.getBaseName().replaceFirst(/info./, "").replaceFirst(/.txt/, ""), it] }
			.set { ch_iPCRinfo }
		
		}
	
		// Map iPCR reads to graph genome with giraffe and sort gam files (will prob need to split to use for counts tables and reads of interest)
		println "Aligning reads"
		ch_unsorted = GIRAFFE(ch_iPCRreads).alignments
		ch_alns = GAM_SORT(ch_unsorted)
		
	}
	else	{
	// Find the info files
		params.infopattern = "${params.storeDir}/${params.outIpcrTrim}/info*txt.gz"
		Channel
			.fromPath(params.infopattern, checkIfExists: true)
			.map { it -> [it.getBaseName().replaceFirst(/info./, "").replaceFirst(/.txt/, ""), it] }
			.set { ch_iPCRinfo }
	
	// Find existing alignments and assumed index locations
		println "Finding aligned reads"
		params.alnpattern = "$params.storeDir/$params.refName/$params.outAln/*.gam"
		Channel
			.fromPath(params.alnpattern, checkIfExists: true)
			.map { it -> [it.getBaseName(), it, "${it}.gai"] }
			.set { ch_alns }
	}
	
	// Option to run giraffe only when space is limited so that gam files can be stashed elsewhere and the rest of the pipeline done as a second run
	if ( !params.giraffeOnly )		{

		// PART 2 - Process gam files and cDNA fastq to get counts tables
		// Get a channel of subgraphs for splitting up the alignment
		params.sgpattern = "$params.subgraphsPath/*vg"
		Channel
			.fromPath(params.sgpattern, checkIfExists: true)
			.map { it -> [it.getBaseName(), it] }
			.set { ch_subgraphs }
		
		if (params.doCountTable)	{
			// Join alignment and info files. If split, this was done before trimming so each aligned fastq should have its own info file
			ch_alns
				.join(ch_iPCRinfo, failOnMismatch: true)
				.set { ch_alninfo }
			
			// Combine subgraph and alignments to get channel with all combinations of subgraph and alignment file
			// Each element should be tuple of chrID, subgraph path, readsID, gam path, gam index path, info path
			ch_subgraphs
				.combine(ch_alninfo)
				.set { ch_tosplit }
			
			// Pass tuple of chr and alignment to process that will extract the gam for the chr and process into iPCR table. Possibly limit number of jobs further than in the eddie config file. Emits tuple of chrID, path to table
			IPCR_GAM(ch_tosplit)
			// Then collect iPCR tables for each chromosome and combine to get channel of chrID-mergediPCRtable
			MERGE_IPCR(IPCR_GAM.out.onetable.groupTuple())
			ch_iPCRtables = MERGE_IPCR.out.tables
	
			// sub-workflow counts
			// For each chromosome combine iPCR table and (genome-wide) cDNA counts (1 per experiment)
			COUNTS(ch_iPCRtables)
			ch_counts_tables = COUNTS.out.counts_table
			
			// Find how many unique fragments and totals of iPCR and cDNA rows
			COUNT_OUTPUT(COUNTS.out.counts_table)
		}
		else	{
			params.ctpattern = "$params.storeDir/$params.refName/$params.outCountsTables/$params.libID*counts_table.txt.gz"
			params.ctreplace = "${params.libID}."
			Channel
				.fromPath(params.ctpattern, checkIfExists: true)
				.map { it -> [it.getBaseName().replaceFirst(/.counts_table.txt/, "").replaceFirst(params.ctreplace,""), it] }
				.set { ch_counts_tables }
		}
		
		// PART 3 - Process counts tables and variants from cactus vcf to find fragments to analyse to find expression differences
		// Get variants and walks, then pair up with counts table batches
		// Assume within the directory the variants are named just with their chromosome IDs .txt
		Channel
			.fromPath("${params.multiDir}/*")
			.map {it -> [it.getBaseName().replaceFirst(/\..*/, ""), it] }
			.set { ch_multi }
		Channel
			.fromPath("${params.varsDir}/*")
			.map {it -> [it.getBaseName().replaceFirst(/\..*/, ""), it] }
			.join(ch_multi, failOnMismatch: true)
			.set { ch_vars }
		
		// Join up with counts table so each has the files of variants for its originating chromosome. If counts table (i.e. chromosome) has no variants file then it won't be in the output channel
		// Identifier is just the counts table name minus file extension
		ch_counts_tables
			.map { it -> [it[0], it[1].getBaseName().replaceFirst(/.txt/,""), it[1]] }
			.combine (ch_vars, by:0)
			.set { ch_counts_vars }
		
		// Find variants in the reads of each counts table (all variants together)
		FIND_VARS(ch_counts_vars)
		
		// Get walks files for imputation to find variants in unsequenced parts of the fragment
		Channel
			.fromPath("${params.walksDir}/*")
			.map {it -> [it.getBaseName().replaceFirst(/\..*/, ""), it] }
			.set { ch_walks }
		
		// Split annotated counts table files into batches of params.fragBatchSize
		// For each batch produces a tuple with chrID, counts table ID, batch filename, two variants files
		// Substitute batch ID for counts table ID then join to walks file
		
		FIND_VARS.out.annot_countstable
			.splitText( by: params.fragBatchSize, file: true, elem: 2, keepHeader: true)
			.map {it -> [it[0], it[2].getBaseName(), it[2], it[3], it[4]] }
			.combine (ch_walks, by:0)
			.set { ch_vtable_walks }
		IMPUTE_VARS(ch_vtable_walks)
		
		// Merge annotated counts tables back into chromosomes and split into one line per variant to give the final output
		ch_varcounts = IMPUTE_VARS.out.imputed_countstable.groupTuple()
		SPLIT_VARS(ch_varcounts)

	}
}
