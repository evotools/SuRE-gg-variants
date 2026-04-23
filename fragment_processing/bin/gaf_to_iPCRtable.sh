#!/usr/bin/env bash

# gaf_to_iPCRtable.sh

# Process gaf file produced by aligning short reads to graph genome with giraffe, converted to gaf and sorted by read name
# Similar to bam_to_iPCRtable.sh but with many edits to parse the different input
# To avoid issues in later script, this version does not have a header line. Instead the header is output as a separate file in case it's ever needed.

# Because the alignment position is given relative to the nodes that make up the path, it's not straightforward to get the strand for the fragment. Originally tried to look at the < and > to guess which read was which but this doesn't work when the assembly contig has been aligned (or assembled, depending on how you look at it) as the reverse complement
# So simply put the F read path as read1 and the R read path as read2
# For strand put NA as I need this column kept to work with everything downstream

# Some read pairs have one read aligned ok and then the other with a very short alignment, to the same node and in the same direction. Exclude such pairs and ones aligned to adjacently-numbered nodes if there are fewer than 20 base matches (column 10 in the GAF). Still checking this although no longer setting F and R

# Only the mapping quality for the first read is given, since it's almost always the same for F and R

# No sequences are included

# inputs are gaf file, info file, name for output

# Ensure that the gam file has been filtered so that lines where only F or R are included have been removed (otherwise will fail at line 93)

GAWK=gawk
CHR="NA"	# If processing from one chromosome can give the ID for the output
OPTIND=1
OPTSTRING="f:i:o:g:c:"
while getopts ${OPTSTRING} opt; do
	case $opt in
	f)
		GAF=${OPTARG}
		;;
	i)
		INFO=${OPTARG}
		;;
	o)
		OUTNAME=${OPTARG}
		;;
	g)
		GAWK=${OPTARG}
		;;
	c)
		CHR=${OPTARG}
		;;
	?)
		echo "Invalid option: -${OPTARG}"
		exit 1
		;;
	esac
done

OUTPUT=${OUTNAME}.iPCR.txt
STATS=${OUTNAME}.stats.txt
HEADER=${OUTNAME}.header
SKIPPED=${OUTNAME}.skipped_pairs.txt

# Make empty file to avoid nextflow error if no pairs are skipped
touch ${SKIPPED}

# For sorting, have an extra column with the first node (then sort on chr and this)

cat ${GAF} | \
${GAWK} -v chromosome=${CHR} -v skipped=${SKIPPED} '
BEGIN { 
  FS  = "\t";
  OFS = "\t";
};
{
	# Check that the first read ends in /1
	if (!($1 ~/\/1$/)) {
		print "Read is not first in pair: "$1  >> "/dev/stderr"
		exit 1
	}
	
	# Check read is correctly formatted
	# If first read of file is incorrectly formatted then exit with error
	# otherwise skip this and the next line to go to the next pair
	if (! ($18 ~/^pd/) ) {
		print "Unexpected gaf formatting: expect pd field in column 18\n"$0 >> "/dev/stderr"
		if (NR == 1)	{
			exit 1
		}
		getline;
		next;
	}
	# Check this read (1) is labelled as being in a proper pair, otherwise skip this and
	# the next read
	if (! ($18 ~/1$/) )	{
		getline;
		next;
	}
	
	# Store alignment data of the first read in the pair
	read[1, "ID"]    = $1;
    read[1, "START"] = $8;
    read[1, "END"] = $9;
    read[1, "PATH"] = $6;
    read[1, "MAPQ"] = $12;
    read[1, "MATCHES"] = $10

	# Read second read in pair
	getline;
	
	# check formatting, skip pair if not ok
	if (! ($18 ~/^pd/) ) {
		print "Unexpected gaf formatting: expect pd field in column 18\n"$0 >> "/dev/stderr"
		next;
	}
	
	read[2, "ID"]    = $1;
    read[2, "START"] = $8;
    read[2, "END"] = $9;
    read[2, "PATH"] = $6;
    read[2, "MAPQ"] = $12;
    read[2, "MATCHES"] = $10
    
    # Check reads actually match
    read1 = read[1, "ID"]
    read2 = read[2, "ID"]
    sub(/\/1$/,"",read1)
    sub(/\/2$/,"",read2)
    if (read1 != read2)	{
    	print "Read IDs do not match: check input has only proper pairs" > "/dev/stderr"
    	print "Read 1 ID: "read[1, "ID"]", Read 2 ID: "read[2, "ID"] > "/dev/stderr"
    	exit 1
    }
    
    # Is fragment (likely to be) on forward or reverse strand relative to the reference?
    # Check for reads where only one really has a decent match (see above)
    if (!(((read[1, "PATH"] ~/>/) && (read[2, "PATH"] ~/</)) || ((read[2, "PATH"] ~/>/) && (read[1, "PATH"] ~/</))))	{
    	# Reads are aligned in the same orientation (should have one with > and one with <)
    	# If also on same node and one is short, exclude; otherwise ok
    	if ((read[1, "MATCHES"] < 20) || (read[2,"MATCHES"] < 20)) {
    		split (read[1,"PATH"], nodes1, "[<>]")
    		split (read[2,"PATH"], nodes2, "[<>]")
    		for (n in nodes1)	{
    			if (n != 1) {
	    			nodes1ID[nodes1[n]] = 1
	    		}
    		}
    		for (n in nodes2)	{
    			if (n != 1)	{
	    			if ((nodes2[n] in nodes1ID) || ((nodes2[n] + 1) in nodes1ID) || ((nodes2[n] - 1) in nodes1ID))	{
    				#	print "Skipping poor pair alignment" > "/dev/stderr"
    					print read[1,"ID"] ","read[1,"PATH"]","read[2,"ID"] ","read[2,"PATH"] >> skipped
    					next;
    				}
    			}
    		}
    	}
    }
    
    # Find the F read to help sorting (may not be the first node in the sequence but just need something here)
    split (read[1, "PATH"], nodes, /<|>/)
    
    # Print the output in the following format (extra column compared with bam version)
    # readID chromosome first.node F.read.pos R.read.pos strand MAPQs
    # MAPQs for F then R reads separated by comma
    print ( read[1, "ID"],
    chromosome,
    nodes[2],
    read[1, "PATH"]"_"read[1, "START"]"-"read[1, "END"],
    read[2, "PATH"]"_"read[2, "START"]"-"read[2, "END"],
    "NA",
    read[1, "MAPQ"])
}
' > ${OUTPUT}

if [ $? != 0 ]; then
  echo "Read parsing failed, check error log"
  exit 1
fi

# Add barcode from INFO file to output file and discard any reads without valid barcodes
# Very similar to the corresponding part of the code in bam_to_iPCRtable.sh
# Then sort on location (chr then node of first read)

mv ${OUTPUT} ${OUTPUT}.tmp

${GAWK} -F '\t' -v statsfile=${STATS} ' 
BEGIN {
	incl=0;
	excl=0;
}
FNR==NR{
	if (FILENAME != ARGV[1]) {exit} # if 1st input file is empty, abort
	a[$1]=$0; # add line to array a with key read ID and value entire line
	next
}
$2==-1 { 
  	# here, the info file indicates the read did not contain an identifiable
  	# barcode. Thus this read should be skipped and removed from the array
  	# Info file is for all reads but only looking at a subset, only count if in array
  	
  	if ($1 in a)	{
		excl++
		delete a[$1]
	}
	next; 
}
{ 
	# Check if the barcode is valid
	if ($1 in a) {
		BC = $5
    	len  = length(BC)
    	BClen[len]++
    	hasN = BC~/N/
    	NNN[hasN]++
    	if (len==20 && !hasN ) {
      		incl++
      		a[$1]=a[$1] FS BC;
	    }
    	else	{
	      excl++
    	  delete a[$1]
		}
	}
}
END {
#  if (FNR == NR) {print "2nd file empty" >> statsfile; exit} # the 2nd input file is empty; abort
  for (key in a) {cnt++; print a[key];}
  print "while filtering for proper barcodes: included = "incl", discarded = "excl >> statsfile
  
  printf "\n\n" >> statsfile
  
  printf ("FragmentCount\t%d\n", cnt) >> statsfile
  
  printf "iPCR_BC_lengths" >> statsfile
  for (k in BClen) printf("\t%d", k) >> statsfile
  printf "\n" >> statsfile
  
  printf "iPCR_BC_lengths_counts" >> statsfile
  for (k in BClen) printf("\t%d", BClen[k]) >> statsfile
  printf "\n" >> statsfile
  
  printf "iPCR_BC_NNN" >> statsfile
  for (k in NNN) printf("\t%d", k) >> statsfile
  printf "\n" >> statsfile
  
  printf "iPCR_BC_NNN_counts" >> statsfile
  for (k in NNN) printf("\t%d", NNN[k]) >> statsfile
  printf "\n" >> statsfile
  
  printf "\n\n" >> statsfile
} ' ${OUTPUT}.tmp <(gzip -cd ${INFO}) | \
# remove read ID column
cut -f 2- |\

# Sort resulting output on chr and first nodeID
sort -k1,1 -k2,2n | \
# remove duplicates and add a count column
uniq -c |\
# reorder columns: move count to 7th column
# and write to $OUTPUT the following columns (1stnode not really needed but placeholder for 'length' in the bam version:
# chr first.read.pos second.read.pos 1stnode strand barcode count MAPQ.Fread
${GAWK} -v statsfile=${STATS} -v headerfile=${HEADER} ' 
BEGIN{ 
  OFS="\t";
  print("chrname","read1_pos","read2_pos","first_node","strand","barcode","count","MAPQs") > headerfile;
} 
{ 
  print($2,$4,$5,$3,$6,$8,$1,$7)
}
END {
  printf("FragmentUniqCount_ReadSequenceIgnored\t%d\n\n", NR) >> statsfile
}' | gzip > ${OUTPUT}.gz

#rm ${OUTPUT}.tmp