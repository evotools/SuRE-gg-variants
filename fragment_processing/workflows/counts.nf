/*
Do steps to go from iPCR table files and raw cDNA reads to counts tables and SuRE score tables
For cDNA reads
	- trim off adapters leaving just the barcodes
	- Filter out invalid barcodes and count, merging together barcodes from multiple sequencing lanes
Then merge with per-chromosome iPCR tables
If cDNA fastq files have already been processed, will need to link to the processed cDNA_counts.txt.gz in the storage directory as this is where the pipeline will look for them
*/

// Import modules
include { CUTADAPT_CDNA } from '../modules/cutadapt.nf'
include { COUNTBARCODE } from '../modules/cDNA_counts.nf'
include { MERGE_IPCR_CDNA; AGG_COUNTS } from '../modules/counts_tables.nf'

workflow COUNTS {
	take:
		ch_iPCRtables
	
	main:
		// Create channel of cDNA libraries (experiments) - for each library, tuple of experiment ID, directory with cDNA fastq files, whether it needs processing to counts (false if counts file already exists)
		// Separate based on whether the cDNAs need processing or not
		Channel.fromPath(params.cDNAList)
			.splitCsv(header:true)
			.map{ row-> tuple("$row.ExpID"),("$row.cDNADir"),("$row.ifCount") }
			.branch	{
				tocount: it[2] == 'true'
				nocount: true
			}
			.set {ch_cDNAs}
	
		// Process cDNA to barcode counts
		// Get read names as channel of [label of ExpID, file path, filename]
		ch_cDNAseq = ch_cDNAs.tocount.map { it -> [it[0], files("${it[1]}/*fq.gz")] }.transpose().map { it -> [ it[0], it[1], it[1].getBaseName().replaceFirst(/.fq/, "") ] }

		// Trim reads (per fastq)
		CUTADAPT_CDNA(ch_cDNAseq)
	
		// Concatenate reads info files and run count_cDNA_barcodes.sh to get counts
		COUNTBARCODE(CUTADAPT_CDNA.out.groupTuple())
		ch_barcodes = COUNTBARCODE.out.counts
		
		// Gather existing counts files and add to newly processed to get channel of expID-cDNAcountsfile pairs
		// Reformat as a string consisting of a space separated string of "ID,file"
		ch_existing = ch_cDNAs.nocount.map { it -> [it[0], file("${params.storeDir}/${params.outCdna}/${it[0]}.cDNA_counts.txt.gz", checkIfExists: true)] }
		ch_barcodes_list = ch_barcodes.concat(ch_existing)
		ch_barcodes_count = ch_barcodes_list.count()
		ch_barcodes_all = ch_barcodes_list	.map { it -> "${it[0]},${it[1]}"}
											.reduce { a, b -> return "$a $b" }
	
		// Combine with iPCR libraries to get a counts table for each
		ch_iPCRtables
			.combine(ch_barcodes_all)
			.set { ch_iPCR_cDNA }
		// Merge iPCR and cDNA counts to get unsorted counts table and library sizes
		ch_counts_all = MERGE_IPCR_CDNA(ch_iPCR_cDNA)
	
		// Sort by position and aggregate positional replicates
		AGG_COUNTS(ch_counts_all.table, ch_barcodes_count)
	
		// Normalise to get SuRE score tables?
	
	// Emit outputs for use in main workflow
	// Channel of chromosome ID-path to counts table
	emit:
		counts_table = AGG_COUNTS.out.counts_table
	//	lib_sizes = AGG_COUNTS.out.lib_sizes
}
