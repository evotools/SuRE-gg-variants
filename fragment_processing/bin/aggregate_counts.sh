#!/usr/bin/env bash

# aggregate_counts.sh

# Original counts table has one line per barcode. This script sorts the counts table by position and aggregates (i.e. sums cDNA counts for) barcodes with identical positions
# Code mostly taken from merge-iPCR-cDNA-plDNA.bash by Ludo Pagie from SuRE-variants-pipeline
# However the 'position-label' has been revised to compare the chr, start, end and read paths of the fragments
# i.e. not looking at where the read ends, and ignoring the strand because this isn't set, and if the nodes match they are on the same strand
# Like with merge_iPCRtables.sh, outputs a list of barcodes that have been merged but were not identical (using the barcode that ends up in the output!)
# In the output use the paths from whichever has the highest iPCR count (or if equal the first in the array is used, which is arbitrary if not random)
# So now will only work with counts tables from the graph genome

# inputs are unsorted counts table, output, sort pattern if not default

# Only the first encountered barcode will be kept. This means when looking up reads by barcode, reads with the subsequent barcode will not be found. However the first read will have the counts for all instances of this DNA segment.

# For iPCR tables from graph genomes, the 5th column will be the ID of the first node of the fragment. If from a linear genome (bam file) it will be the fragment length. Want to include this for sorting for graph genomes but not for linear. So if iPCR table came from a bam file, use the following as the sort pattern (set as nextflow parameter countTableSort in in pipeline):
# "-k2,2 -k3,3n -k4,4n"
# If using chromosomes with a fixed prefix (e.g. 'chr1') can use (e.g.) -k2.4,2 to start at the 4th character of the 2nd field to save (possibly minimal) time
# The sort pattern used below for the graph genomes is not perfect as it sorts the paths alphabetically, but it will be enough to get the rows together that need to be together

GAWK=gawk
SORTPATTERN="-k2,2 -k5,5 -k3,3 -k4,4 -k6,6" # Default for alignments from graph genome by chromosome, first node, F read path, R read path, strand

OPTIND=1
OPTSTRING="u:o:s:g:n:"
while getopts ${OPTSTRING} opt; do
	case $opt in
	u)
		UNSORTED=${OPTARG}
		;;
	s)
		SORTPATTERN=${OPTARG}
		;;
	o)
		OUTNAME=${OPTARG}
		;;
	g)
		GAWK=${OPTARG}
		;;
	n)
		NSAMPLES=${OPTARG}
		;;
	?)
		echo "Invalid option: -${OPTARG}"
		exit 1
		;;
	esac
done

# create empty pos_multi file to avoid an error if there are no reads to write to it
touch ${OUTNAME}_pos_multi_BC.txt
touch ${OUTNAME}_merged.txt
touch ${OUTNAME}_merged_nonidentical.details.txt

gzip -cd ${UNSORTED} | \
${GAWK} -v SORTP="${SORTPATTERN}" ' BEGIN {OFS="\t"} 
	NR==1	{
		print $0;
		next
	} # keep header
	{
		print $0 | "sort "SORTP
	}' | \
${GAWK} -v POS_MULTI_BC_FNAME="${OUTNAME}_pos_multi_BC.txt" -v NUM_SAMPLES="${NSAMPLES}" -v MULTI_NOTID="${OUTNAME}_merged.txt" -v MULTI_NOTID2="${OUTNAME}_merged_nonidentical.details.txt" '
# awk script to merge counts for duplicated positions
# for every input line a position-label is created as a concatenation of the fields chr/start/end/strand
# if position-labels of subsequent lines are identical the lines are merged;
# ie. the counts for all samples in all lines are summed
# Barcode of the first line is retained

function PROC_PREV_POS(position) {
  # function to process the previous collected lines which have the same position-label
  	if (length(PREVLINE) == 1) {
		print PREVLINE[1]
    	return
 	}
 	# Get output line from first line details and sum counts columns
	split(PREVLINE[1], outline)
	# set counts for iPCR and sample columns to 0
	for (i=8; i<=8+NUM_SAMPLES; i++)
	  outline[i] = 0
	# collect MAPQs to merge
	mapq_str = outline[7]
	split("",mapq)
	split(mapq_str, mapq1, ",")
	for (i in mapq1)	{
		mapq[mapq1[i]] = 1
	}
	# collect complete paths and possible output columns for log
	split("",readpaths)
	split("",pathbarcodes)
	split("",ifkept)
	for (line in PREVLINE) {
	  print PREVLINE[line] > POS_MULTI_BC_FNAME
	  split(PREVLINE[line], w)
	  ifkept[line] = w[1]"\t"w[3]"\t"w[4]"\t"w[7]"\t"w[8]
	  split(w[7], rowmapq, ",")
	  for (one in rowmapq) {
  	    if (!(rowmapq[one] in mapq)) {
		  mapq[rowmapq[one]] = 1
		  outline[7] = outline[7] "," rowmapq[one]
		}
	  }
	  for (i=8; i<=8+NUM_SAMPLES; i++)
		outline[i] = outline[i] + w[i]
	  fullpos = w[3]"\t"w[4]
	  if (fullpos in readpaths)	{
	  	readpaths[fullpos] = readpaths[fullpos] + w[8]
	  }
	  else	{
	  	readpaths[fullpos] = w[8]
	  	pathbarcodes[fullpos] = w[1]
	  }
	}
	# if alignments are not identical then choose whichever has the highest iPCR count
	# output details of all to log file
	if (length(readpaths) > 1)	{
	  readpos_out = ""
	  maxcov = 0
	  for (pos in readpaths) {
	    if (readpaths[pos] > maxcov) {
	      maxcov = readpaths[pos]
	      readpos_out = pos
		}
	  }
	  split(readpos_out,pos_array)
	  outline[3] = pos_array[1]
	  outline[4] = pos_array[2]
	  outline[1] = pathbarcodes[readpos_out]
	  
	  # log
	  for (det in ifkept) {
	    print ifkept[det]"\t"outline[1]"\t"outline[3]"\t"outline[4] > MULTI_NOTID2
	  }
	}
	
	# print the merged output line
	printf("%s", outline[1])
	for (i=2; i<=length(outline); i++)
	  printf("\t%s", outline[i])
	printf("\n")
	
	# print the output barcode and number of merged barcodes to another log file
	print outline[1]"\t"length(PREVLINE) > MULTI_NOTID
		
}

BEGIN {
  FS="\t"
  OFS="\t"
}

NR == 1 { # read and print header from input
  print
}

NR == 2 { # initialize PREVPOS and PREVLINE with values from 1st input line
  POS1 = $3
  gsub(/-[0-9]+$/, "", POS1)
  POS2 = $4
  gsub(/-[0-9]+$/, "", POS2)
  PREVPOS = $2"\t"POS1"\t"POS2
  PREVLINE[1] = $0
}

NR > 2 {
  POS1 = $3
  gsub(/-[0-9]+$/, "", POS1)
  POS2 = $4
  gsub(/-[0-9]+$/, "", POS2)
  POS = $2"\t"POS1"\t"POS2
  if (POS == PREVPOS)	{
    PREVLINE[length(PREVLINE)+1] = $0
  }
  else {
    PROC_PREV_POS(PREVPOS)
    delete PREVLINE
    PREVLINE[1] = $0
    PREVPOS     = POS
  }
}

END {
  # for sensible error message if no reads
  if (length (PREVLINE) == 0 )	{
  	print "No read lines"
  	exit 1
  }
  PROC_PREV_POS()
} ' | gzip -c > ${OUTNAME}.gz

if [ ${PIPESTATUS[2]} != 0 ]; then
  echo "Likely no reads aligned for $UNSORTED"
  exit 1
fi


gzip ${OUTNAME}_pos_multi_BC.txt
gzip ${OUTNAME}_merged*.txt
