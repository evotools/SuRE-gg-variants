#!/usr/bin/env bash

# merge_iPCRtables.sh
# Merge the output of aln-to-iPCRtable from multiple sequencing lanes for the same sample
# Count duplicates (same barcode, same mapping) and reduce to one
# Update: for barcodes where the read alignment paths do not match exactly, instead of excluding completely as in the original version, keep if the only difference is the exact end position of the read alignment, and the read start positions and the nodes passed through are identical. In this case take the path and other details of the read with the highest MAPQs

# Mostly taken from iPCR-merge-bedpe-Filter-BC-multi-pos.bash by Ludo Pagie from SuRE-variants-pipeline

# Works on iPCR tables generated from either gam/gaf alignments, after adding code to parse the read paths this no longer works for iPCR tables from bam files

# Pass input filenames to script as arguments

GAWK=gawk

OPTIND=1
OPTSTRING="o:c:g:"
while getopts ${OPTSTRING} opt; do
	case $opt in
	o)
		OUTNAME=${OPTARG}
		;;
	c)
		NCORES=${OPTARG}
		;;
	g)
		GAWK=${OPTARG}
		;;
	?)
		echo "Invalid option: -${OPTARG}"
		exit 1
		;;
	esac
done

shift $(( OPTIND - 1 ))

BEDPE_FNAME=$*

touch merge-bedpe-log.txt
touch multi-end-pos.txt

# Prints barcodes with multiple lines to log file for checking. Because this is in the nextflow pipeline I'm just using a fixed filename that will be in the workdir
gzip -cd ${BEDPE_FNAME} | \
  # sort bedpe data on barcode (dictionary order)
  sort --parallel=${NCORES} -k 6,6d | \
  ${GAWK} -v BC_MULTI_POS_FNAME="merge-bedpe-log.txt" -v BC_MULTI_KEPT="multi-end-pos.txt" '
  # this piece of awk collapses iPCR bedpe data "points": 
  # while consecutive lines contain identical barcodes the lines are collected in a array "PREVLINE"
  # once a line with a different barcode is read the "PREVLINE" array is
  #   processed before resetting the array with the new line as only element
  # processing of the array means that all bedpe data "points" are grouped on
  #   genomic position; the position with max readcount is chosen as "true"
  #   barcode; this barcode is printed with corrected values for read-count and
  #   mapq
function GETSTART(path) {
	loc=""
	match(path, /[<>][0-9]+/, nodes)
	match(path, /_([0-9]+)-/, pos)
	loc=nodes[0]"_"pos[1]
	return loc
}
function PROC_PREV_BC()	{
  	# function to process previously collected barcodes (to avoid a bug that missed the last line)
  	if (length(PREVLINE) == 1) {
		# if previous barcode appears just once then just print the corresponding line
    	print PREVLINE[1]
	} else {
		# if the barcode appears multiple times then may be mapped to multiple positions so check
		# Skip unless fragment start and end and the nodes passed through are the same in all reads
		# If read ends differ, use result for read with best MAPQ (choose at random if >1)
		# If kept, all reads count towards iPCR total
		delete BC2POS
		count=0
		split("",mapqs)
		for (line in PREVLINE) {
			split(PREVLINE[line], records)
			count += records[COLCNT]
			POS = records[COLCHR]"\t"records[COLS]"\t"records[COLE]"\t"records[COLLENGTH]"\t"records[COLSTR]
			if (!(POS in BC2POS))	{
				BC2POS[POS] = BC2POS[POS]
			}
			mapqs[line]=records[COLMAPQ]
			print PREVLINE[line] > BC_MULTI_POS_FNAME
		}
		# if > position string then check if can still keep it
		bestline=""
		if (length(BC2POS) > 1)	{
			delete FRAGPOS
			delete BCNODES
			for (line in PREVLINE) {
				split(PREVLINE[line],cols)
				FRAGSTART=GETSTART(cols[2])
				FRAGEND=GETSTART(cols[3])
				FRAGBOTH=FRAGSTART FRAGEND
				if (!(FRAGBOTH in FRAGPOS))	{
					FRAGPOS[FRAGBOTH] = FRAGBOTH
				}
				split(cols[2],startnodes,"_")
				split(cols[3],endnodes,"_")
				NODESBOTH=startnodes[1] endnodes[1]
				if (!(NODESBOTH in BCNODES))	{
					BCNODES[NODESBOTH] = NODESBOTH
				}
			}
			if (length(FRAGPOS) > 1 || (length(BCNODES) > 1))	{
				print "###\n" > BC_MULTI_POS_FNAME
				return
			}
			
			# Choose a line with highest mapq
			bestline=""
			maxmapq=(-50)
			for (line in mapqs) {
				if (mapqs[line] > maxmapq)	{
					maxmapq=mapqs[line]
				}
			}
			split("",best)
			for (line in mapqs) {
				if (mapqs[line] == maxmapq)	{
					best[length(best)+1]=line
				}
			}
			bestline=best[ 1 + int(rand()*length(best)) ]
			print cols[COLBC]"\t"length(BC2POS) > BC_MULTI_KEPT
		}
		else	{
			bestline=1
		}
		
		# Check MAPQs and output all unique (for filtering later)
		mapq_str = ""
		split("",mapq_seen)
		for (line in mapqs)	{
			if (!mapq_seen[mapqs[line]])	{
				mapq_seen[mapqs[line]] = 1
				if (mapq_str) {
					mapq_str = mapq_str "," mapqs[line]
				}
				else	{
					mapq_str = mapqs[line]
				}
			}
		}
		
		# Use either the first line or the best MAPQ line (if alignments different)
		# Update counts to sum all lines, and MAPQ
		split(PREVLINE[bestline], output)
		output[COLCNT] = count
		output[COLMAPQ] = mapq_str
		# print output
		outline=output[1]
		for (i=2; i<=length(output); i++)	{
			outline=outline "\t" output[i]
		}
		print outline
		# and print to logfile with # to indicate it was printed to output
		print "#"outline"\n" > BC_MULTI_POS_FNAME
		
	}
}
  
BEGIN {
    FS="\t"
    OFS="\t"

    # set some column numbers
    COLLENGTH=4
    COLMAPQ=8
    COLBC=6
    COLCHR=1
    COLS = 2
    COLE = 3
    COLSTR = 5
    COLCNT = 7
}
NR == 1 {
	PREVBC      = $COLBC
    PREVLINE[1] = $0
    next
}
{
	# here the fragment is considered good quality (ie properly aligned and such).
    # the fragments are sorted on barcode. If consecutive barcodes are
    # identical all the counts are summed and one line is output
#	print "processing line " NR" "$COLBC" "PREVBC
	if ($COLBC == PREVBC) {
      # we read an identical barcode; add this line to the line stack
      PREVLINE[length(PREVLINE)+1] = $0
    } else {
		# we have come to a novel barcode; 
		# process the previous barcode
    	PROC_PREV_BC()
    	
    	# previous barcode is dealt with. start a new line stack
		delete PREVLINE
		PREVLINE[1] = $0
		PREVBC = $COLBC
    }
}
END	{
	PROC_PREV_BC()
}  ' | gzip > ${OUTNAME}.gz

