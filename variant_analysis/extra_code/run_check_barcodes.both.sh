#!/bin/sh
#$ -N mmapBoth
#$ -cwd
#$ -l h_vmem=24G
#$ -l h_rt=24:00:00
#$ -P roslin_prendergast_cores
#$ -pe sharedmem 16
#$ -R y

# Similar to run_check_barcodes but looks at both at the same time
# Link all counts tables to wd

# Directory containing the *.counts_table.txt.gz files
DIR=$1

# Set number of threads to 8
THREADS=16

cd "$DIR" || exit 1

# If you have pigz installed, uncomment this for multi-threaded decompression:
DECOMP="pigz -dc"
# Otherwise, use gzip:
#DECOMP="gzip -dc"

# 1) Decompress all files, extract barcode (col 1) and chromosome (col 2)
# 2) Sort by barcode using multi-threaded GNU sort
# 3) Awk pass over sorted data to build per-barcode chromosome list
#    and print only those barcodes seen more than once
# 4) Second awk pass over duplicated barcodes to print each to appropriate per-chromosome file

find . -name '*counts_table.txt.gz' -print0 \
  | xargs -0 $DECOMP \
  | awk -F'\t' '{print $1"\t"$2}' \
  | sort --parallel="$THREADS" -S 80% -k1,1 \
  | awk -F'\t' '
    {
      total++                           # count every line (all copies, including duplicates)

      if ($1 == prev_barcode) {
        chroms = chroms "," $2
        count++
      } else {
        # New barcode: emit previous if it was duplicated
        if (NR > 1 && count > 1) {
          print prev_barcode "\t" chroms "\t" count
        }
        prev_barcode = $1
        chroms = $2
        count = 1
      }
    }
    END {
      # Handle last barcode
      if (count > 1) {
        print prev_barcode "\t" chroms "\t" count
      }
      # Output total number of lines (including all copies of duplicates)
      print "#TOTAL_LINES\t" total > "total_line_count.txt"
    }
  ' > "bothlibs.duplicated_barcodes_with_chroms.tsv"

# Second pass to get by-chromosome lists
gawk -F'\t' '
{
	split($2, chroms, ",")
	for (i in chroms) {
    	print $1 >> "both." chroms[i] ".dup_bc.tsv"
  	}
  	
}' "bothlibs.duplicated_barcodes_with_chroms.tsv"