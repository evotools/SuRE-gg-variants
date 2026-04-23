#!/usr/bin/env bash

# gam2bam_to_iPCRtable.sh
# Version of bam_to_iPCRtable.sh that's been tweaked to be able to work with the bam files produced by surjecting the gam files from the main pipeline.
# Also updated to avoid the bug in line 160 when using just chromosome numbers as the ID

# Main difference in the bam file from surject is that the QC flags are missing so all will fail the checks for proper pairs and first in pair.
# Pre-filter for proper pairs when extracting from the gam file. This is not as stringent as bowtie would use but I think it will do
# To check for the appropriate read being the first in the pair, just check it ends in /1 (as in the gaf version of the script)
# The flag for read on R strand is present ()

# Process bam file of short reads aligned to linear genome to make the OUTPUT-like file used to generate the count tables
# Mostly taken from iPCR-bam2OUTPUT.bash and iPCR-map-BC.bash by Ludo Pagie from SuRE-variants-pipeline
# To avoid issues in later script, this version does not have a header line. Instead the header is output as a separate file in case it's ever needed.

# inputs are bam file, info file, name for output, chr name (since bam file has path name and we want chr name)

GAWK=gawk
OPTIND=1
OPTSTRING="b:i:o:g:c:"
while getopts ${OPTSTRING} opt; do
	case $opt in
	b)
		BAM=${OPTARG}
		;;
	i)
		INFO=${OPTARG}
		;;
	o)
		OUTNAME=${OPTARG}
		;;
	c)
		CHR=${OPTARG}
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

OUTPUT=${OUTNAME}.iPCR.txt
STATS=${OUTNAME}.stats
HEADER=${OUTNAME}.header
SKIPPED=${OUTNAME}.skipped_pairs.txt

# Make empty file to avoid nextflow error if no pairs are skipped
touch ${SKIPPED}

samtools view ${BAM} | \
  ${GAWK} -v skipped=${SKIPPED} -v chrid=${CHR} '

function CIGAR2length( cigar ) {
  # compute from a CIGAR string the length of the aligned genomic region:
  # first; CIGAR operators which appear only once do not have a counter in
  # front of them. First I will place a '1' in front of those so afterwards I
  # can easily split the CIGAR string and have counts and operators
  # separately.
  # if first character is not a count but a CIGAR operator (MDINSHP=X), put a 1 in front of it
  mod=gensub("^([MDINSHP=X])", "1\\1",1,cigar);
  # for all consecutive pairs of CIGAR operators insert a 1 as count for the 2nd operator
  mod=gensub("([MDINSHP=X])([MDINSHP=X])", "\\11\\2",1,cigar);
  while (mod != cigar) {
    cigar=mod;
    mod=gensub("([MDINSHP=X])([MDINSHP=X])", "\\11\\2",1,cigar);
  }
  
  # split CIGAR string into an array of the counts and an array of the CIGAR operators
  split(cigar, counts, /[MDINSHP=X]+/);
  split(cigar, ops, /[[:digit:]]+/);
  # map the operators to 0,1; depending on whether the operators corresponds to
  # an increase in length of the aligned genomic region
  for (i in ops) {
    gsub("[MDNPX=]", "1", ops[i]) # these ops increase genomic alignment
    gsub("[ISH]", "0", ops[i]) # these ops do not increase genomic alignment
  }
  # compute the total length by multiplying the counts and the lengths of the operators
  len=0;
  for (i=1; i<length(counts); i++) {
    len = len + (counts[i]*ops[i+1])
  }
  return len;
}

BEGIN { 
  FS  = "\t";
  OFS = "\t";
  flag_readUnmapped = lshift(1, 2);
  flag_readReverseStrand = lshift(1, 4);
};

{
  	# Check that the first read ends in /1
	if (!($1 ~/\/1$/)) {
		print "Read is not first in pair: "$1
	}

    # store alignment data of first read in pair
    read[1, "MAPQ"]  = $5;
    read[1, "START"] = $4;
    read[1, "FLAG"]  = $2;
    read[1, "CIGAR"] = $6;
    read[1, "RNAME"] = $3;
    read[1, "ID"]    = $1;
    read[1, "STRAND"] = "+"
    # compute end position on genome
    read[1, "END"] = read[1,"START"] - 1 + CIGAR2length(read[1, "CIGAR"])

    # read second read in pair and parse data
    getline;
    read[2, "MAPQ"]  = $5;
    read[2, "START"] = $4;
    read[2, "FLAG"]  = $2;
    read[2, "CIGAR"] = $6;
    read[2, "ID"]    = $1;
    read[2, "STRAND"] = "+"
    
    # compute end position on genome
    read[2, "END"] = read[2,"START"] - 1 + CIGAR2length(read[2, "CIGAR"])
    
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
    
    
    # Check that both reads are mapped (surject from gam might not be able to find mapping)
    if ( ( and (read[1, "FLAG"], flag_readUnmapped) ) || ( and (read[2, "FLAG"], flag_readUnmapped) ) )	{
    	print "Read_unmapped: "read[1,"ID"] ","read[1,"FLAG"]","read[2,"ID"] ","read[2,"FLAG"] >> skipped
    	next;
    }
	
	# gam-to-bam conversion even with vg filtering for proper pairs can leave read pairs mapped in the same orientation
	# Check for this and exclude
    if ( and(read[1, "FLAG"], flag_readReverseStrand) ) {
    	read[1, "STRAND"]="-";
    }
  	if ( and(read[2, "FLAG"], flag_readReverseStrand) ) {
    	read[2, "STRAND"]="-";
    }
    if (read[1, "FLAG"] == read[2, "FLAG"])	{
    	print "Same_orientation: "read[1,"ID"] ","read[1,"STRAND"]","read[2,"ID"] ","read[2,"STRAND"] >> skipped
    	next;
    } else	{
    	strand = read[1, "STRAND"]
    }
    
    # But even if in opposite orientations they might be <-- --> instead of proper --> <--, so check this too while working out which is the first_read. Allow overlap but R2 must not start before R1. N.B. Start always < end for the same fragment
    if (strand == "+")	{
    	if (read[1, "START"] > read[2, "START"])	{
    		print "Improper_pos: "read[1,"ID"] "," read[1,"START"] "," read[1,"END"] "," read[1,"STRAND"] "," read[2,"START"] "," read[2,"END"] "," read[2,"STRAND"] >> skipped
    		next;
    	} else {
    		first_read=1
    	}
    } else {
    	if (read[1, "START"] < read[2, "START"])	{
    		print "Improper_pos: "read[1,"ID"] "," read[1,"START"] "," read[1,"END"] "," read[1,"STRAND"] "," read[2,"START"] "," read[2,"END"] "," read[2,"STRAND"] >> skipped
    		next;
    	} else {
    		first_read=2
    	}
    	
    }
    last_read=3-first_read

	# print the output, including the extra data, in the following format:
	# readID seqname start end strand MAPQ.1
	print( read[1, "ID"],
	chrid,
	read[first_read, "START"],
	read[last_read, "END"],
	strand,
	read[first_read, "MAPQ"] )
}
' > ${OUTPUT}

## add barcode sequence from INFO file to OUTPUT file
# $INFO is the info-file created while trimming the forward read by
# cutadapt (ie very first step)
mv ${OUTPUT} ${OUTPUT}.tmp
# merge barcodes from $INFO into $OUTPUT using awk
# (store readID from column-5 from $OUTPUT as key in array 'a' with
# entire line, except readID, as value. read $INFO, look for key and add barcode
# sequence $5) to array element
# info file may contain reads for which no adapter sequence was found: discard
# those reads from the info file prior to adding the barcodes into the OUTPUT
# file

${GAWK} -F '\t' -v statsfile=${STATS} ' 
BEGIN {
  incl=0;
  excl=0;
}
FNR==NR{
  if (FILENAME != ARGV[1]) {exit} # if 1st input file is empty, abort
  a[$1]=$0; # add entire line to array a with key read ID
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
    else {
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

# sort resulting output file on seqname and start (latter numeric
# sort)
sort -k1,1 -k2,2n | \
# remove duplicates and add a count column
uniq -c |\
# reorder columns: move count to 7th column
# and write to $OUTPUT the following columns:
# chr start end length strand barcode count MAPQ
${GAWK} -v statsfile=${STATS} -v headerfile=${HEADER} ' 
BEGIN{ 
  OFS="\t";
  print("seqname","start","end","length","strand","barcode","count","MAPQ") > headerfile;
} 
{ 
  print($2,$3,$4,$4-$3+1,$5,$7,$1,$6)
}
END {
  printf("FragmentUniqCount_ReadSequenceIgnored\t%d\n\n", NR) >> statsfile
}' | gzip > ${OUTPUT}.gz

#rm ${OUTPUT}.tmp