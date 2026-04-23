#!/usr/bin/env bash

# count-cDNA-barcodes.sh
# Parse the info file to get the barcodes, filter out invalid barcodes and count occurrences of each
# Also create files to get stats on barcode length and number of Ns

# Code mostly taken from cDNA-plDNA-count-BC by Ludo Pagie from SuRE-variants-pipeline

# Give all sequencing read files for the experiment on the command line, read all together and sort, so the barcode counts are for the whole experiment (1 replicate of 1 library in 1 cell line) not for the individual sequencing files

## NOT TESTED but used previously in a very similar form

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

INFO_FNAME=$*

NNNFNAME=NNN.${OUTNAME}
BCLENGTHFNAME=BClengths.${OUTNAME}
BCCOUNTFNAME=BCcounts.${OUTNAME}

gzip -cd ${INFO_FNAME} | cut -f 5 | sort --parallel=${NCORES} | uniq -c | \
  gawk -v NNNfname=${NNNFNAME} -v BClengthfname=${BCLENGTHFNAME} '
  BEGIN { 
    OFS = "\t"; 
  } 
  { 
    l = length($2) # length of the barcode sequence
    N = $2 ~ /N/   # boolean specifying whether barcode sequence contains Ns
    if(l==20 && N==0) print $2, $1 # only print correct barcodes
    len[l]++  # collect barcode length statistics 
    NNN[N]++  # collect barcode containing Ns statistics
  }
  END {
    # print the collected statistics of the barcodes to file
    PROCINFO["sorted_in"] = "@ind_num_asc"
    for (i in len) 
      printf("%d\t%d\n", i, len[i]) > BClengthfname
      for (i in NNN)
        printf("%d\t%d\n", i, NNN[i]) > NNNfname
  } ' | \
gzip > ${OUTNAME}.gz

gzip -cd ${OUTNAME}.gz | cut -f 2 | sort | uniq -c | sort -k2g > ${BCCOUNTFNAME}