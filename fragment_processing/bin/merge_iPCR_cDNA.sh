#!/usr/bin/env bash

# merge_iPCR_cDNA.sh
# Merge the files created from the iPCR reads with the barcode counts from one sample
# Give iPCR file as -i
# cDNA barcodes counts files are the remaining options, given as SAMPLEID,FILEPATH (no space)

# Code mostly taken from merge-iPCR-cDNA-plDNA.bash by Ludo Pagie from SuRE-variants-pipeline

# At this stage, identical positions with different barcodes are not aggregated (difference from published). This can happen downstream in a separate script (because different inputs might need sorting differently)

# Output file is different from original
#	- barcode retained as first column, with header (cut -f 2- if want to remove)
#	- Output is not sorted by genome position
#	- Lines with identical position but different barcodes have not been merged

# Also counts the total number of reads in each library in the final counts table. As originally written, this will be identical to the counts in the output of the next script, aggregate_counts.sh, but if any change is made to remove any barcodes from the library at a later stage than this (e.g. filtering for mapping quality) then these library sizes won't be accurate.

# For iPCR tables from graph genomes, the 4th column will be the ID of the first node of the fragment. If from a linear genome (bam file) it will be the fragment length. In either case it will be the 5th column of the output.

GAWK=gawk

OPTIND=1
OPTSTRING="o:c:g:i:"
while getopts ${OPTSTRING} opt; do
	case $opt in
	o)
		OUTNAME=${OPTARG}
		;;
	c)
		NCORES=${OPTARG}
		;;
	i)
		IPCR_FNAME=${OPTARG}
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

TOTALS="lib_sizes.${OUTNAME}"

# From the end arguments, get array of samples and associative array of counts file paths
declare -A FNAMES
declare -a SAMPLES
for pair in $*; do
	IFS=',' read -ra pairing <<< "$pair"
	sample=${pairing[0]}
	fname=${pairing[1]}
	SAMPLES+=("${sample}")
	FNAMES["${sample}"]="${fname}"
done

#### GAWK SCRIPT to combine the files together
{ echo -e "iPCR\t${IPCR_FNAME}";
  for (( i=0; i<${#SAMPLES[@]}; i++ )); do
    echo -e "${SAMPLES[$i]}\t${FNAMES["${SAMPLES[$i]}"]}"; 
  done } | \
${GAWK}  -v TOTALSOUT="${TOTALS}" ' 
BEGIN {
  # global variables
  FS="\t"

  split("", pipes)
  split("", samples)
  split("", totals) # running total of the size of each library

}
NR==1 {
  # iPCR filename
  ipcr_pipe  = "gzip -dc "$2
  ipcr_total = 0
  next
}
{
  # cDNA/plDNA filename
  samples[length(samples)+1] = $1
  pipes[$1]                  = "gzip -dc "$2
  totals[$1] = 0
}
END {
  # print header to stdout 
  printf ("barcode\tchr\tstart\tend\tfirst_node\tstrand\tMAPQs\tiPCR")
  for (i=1; i<=length(samples); i++)
    printf("\t%s", samples[i])
  printf("\n")

  # read input from iPCR datafile; for every BC read from iPCR iterate over all
  # other input pipes and for each pipe extract all lines until the BC read
  # from that input pipe is (alphabetically) larger than the iPCR barcode
  cnt=0
  while ((ipcr_pipe | getline) > 0) {
    # initialize/empty an array to collect fields for output
    split("", lineout)
    BCipcr = $6
    lineout["iPCR"] = $7
    lineout["CHR"]  = $1
    lineout["START"] = $2
    lineout["END"] = $3
    lineout["FIRSTNODE"] = $4	# Or length if linear reference. Keep in for sorting.
    lineout["STRAND"] = $5
    lineout["MAPQ"] = $8
    ipcr_total += $7
    
    # iterate over all open pipes
    for (i in samples) {
      # initialize count for current sample to zero
      lineout[samples[i]] = 0
      # check whether pipe to current sample file is closed; skip to next sample in that case
      if ( !(samples[i] in pipes) ) {
        continue
	  }
	  
	  # for current iPCR barcode check whether the previously read barcode for current sample is identical or not
      if (BCprev[i] == BCipcr)
        # identical barcodes; set count for current iPCR barcode
        lineout[samples[i]] = CNTprev[i]

      if (BCprev[i] > BCipcr) {
        # previous sample barcode is "larger" than iPCR barcode, therefore
        # current sample does not contain current iPCR barcode; continue to
        # next sample
        continue
      }

      # at this point the previous sample barcode is <= current iPCR barcode;
      # we need to read from the current sample pipe until the sample barcode
      # is > iPCR barcode. If a newly read sample barcode is identical to iPCR
      # barcode update the count for the current sample
      while(1) {
        # read from current sample pipe
        # print "in while loop: i="i", sample="samples[i]", pipe="pipes[samples[i]] > "/dev/stderr"
        status = (pipes[samples[i]]) | getline line
        # checkfile read status; if EOF close this pipe and delete the pipe from array _pipes_
        if (status == 0) {
          # print "in while loop, status==0: i="i", sample="samples[i]", pipe="pipes[samples[i]] > "/dev/stderr"
          close (pipes[samples[i]])
          delete pipes[samples[i]]
          # break from _while(1)_ loop which reads from current sample pipe, to read next sample pipe
          break
        }

        # process record read from sample pipe; compare sample barcode and iPCR barcode
        # N.B. My US script swaps the barcode and count fields around cf the original
        split(line, fields)
        if (fields[1] == BCipcr) {
          # sample barcode equal to iPCR barcode; update count for current sample
          lineout[samples[i]] = lineout[samples[i]] + fields[2]
        }
        if (fields[1] > BCipcr) {
          # newly read sample barcode > iPCR barcode, therefore store sample barcode for checking with subsequent iPCR barcodes
          BCprev[i] = fields[1]
          CNTprev[i] = fields[2]
          # break from _while(1)_ loop which reads from current sample pipe, to read next sample pipe
          break
        }
      } # end loop _while(1)_
    } # end loop _for (i in samples)_
    
    # all sample pipes for current iPCR barcode have been processed
    # print record to stdout
    # Format depends on whether iPCR is from linear or graph genome
    if (lineout["START"] ~/[<>]/ || lineout["END"] ~/[<>]/)	{
   	    printf("%s\t%s\t%s\t%s\t%d\t%s\t%s\t%d", BCipcr, lineout["CHR"], lineout["START"], lineout["END"], lineout["FIRSTNODE"], lineout["STRAND"], lineout["MAPQ"], lineout["iPCR"])
	}
	else	{
	    printf("%s\t%s\t%d\t%d\t%d\t%s\t%s\t%d", BCipcr, lineout["CHR"], lineout["START"], lineout["END"], lineout["FIRSTNODE"], lineout["STRAND"], lineout["MAPQ"], lineout["iPCR"])
    }
    
    for (i=1; i<=length(samples); i++)
      printf("\t%d", lineout[samples[i]])
      printf("\n")
    
    # add total to samples total
    for (i=1; i<=length(samples); i++) {
      totals[samples[i]] = totals[samples[i]] + lineout[samples[i]]
    #  print BCipcr" "samples[i]" "totals[samples[i]] >> TOTALSOUT
    }
    
  } # end loop _while ((ipcrpipe | getline) > 0)_
  close (ipcr file)
  # close any pipe left open
  for (p in pipes) {
    # print "pipe "p" still open" > "/dev/stderr"
    # (pipes[p]) | getline > "/dev/stderr"
    close (pipes[p])
    delete pipes[p]
    # print
  }
  # print the totals with headers
  printf("iPCR") > TOTALSOUT
  for (i=1; i<=length(samples); i++)
    printf("\t%s", samples[i]) >> TOTALSOUT
  printf("\n") >> TOTALSOUT
  printf("%d", ipcr_total) >> TOTALSOUT
  for (i=1; i<=length(samples); i++)
    printf("\t%d", totals[samples[i]]) >> TOTALSOUT
  	  
}' | gzip > ${OUTNAME}.gz

