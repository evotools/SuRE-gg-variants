#!/usr/bin/env bash

# project_to_ref.sh
# Given the alignments of the F and R reads of the fragment relative to the graph, find the 
# fragment start and end relative to the reference (simple-minded surject)

# To get the fragment start and end, I need to process the F and R reads in the same way, since
# start of the R read is the end of the fragment

# arguments are
# l - list of barcodes to include (1 per line)
# b - table of barcodes with their start and end alignments from gam
# g - path to subgraph.vg of this chromosome
# v - path to vg binary

# Check in output file that barcode is col1

# Outputs are a table of the barcode, original start and end, and extra columns for the ref start and end positions, and
# a list of errors -- barcodes where either no nodes were found in the reference (M for missing), or where the
# first node wasn't in the reference path (but later ones were) (S for skipped)

OPTIND=1
OPTSTRING="l:b:g:v:p:"
while getopts ${OPTSTRING} opt; do
	case $opt in
	l)
		BCLIST=${OPTARG}
		;;
	b)
		BARCODES=${OPTARG}
		;;
	g)
		GRAPH=${OPTARG}
		;;
	v)
		VGPATH=${OPTARG}
		;;
	p)
		PREFIX=${OPTARG}
		;;
	?)
		echo "Invalid option: -${OPTARG}"
		exit 1
		;;
	esac
done


# Get the relevant barcodes
gawk 'NR==FNR{a[$1]; next} $1 in a' ${BCLIST}  <(gzip -cd ${BARCODES}) > ${PREFIX}.bc_to_check.txt

# Extract the nodes from the alignments and output as node range
gawk '
BEGIN {
	OFS = "\t"
}
function list_nodes(path) {
	delete split_path
	delete nodes
	split(path, split_path, /_/)
	split(split_path[1], nodes, /[<>]/)
	
	nodeslist=""
	for (i = 2; i <= length(nodes); i++) {
		nodeslist = nodeslist nodes[i] "\n"
	}
	return nodeslist
}
{
	printf "%s%s", list_nodes($2), list_nodes($3)
}
' ${PREFIX}.bc_to_check.txt \
 | sort -u > ${PREFIX}.nodes_to_check.txt

mkdir tmp_${PREFIX}

# Extract the nodes from the graph
# If in the reference path, extract the walk lines
$VGPATH chunk -c 0 -R ${PREFIX}.nodes_to_check.txt -O gfa -x ${GRAPH} -b tmp_${PREFIX}/chunk
rp=`$VGPATH paths -L -x ${GRAPH}`
gawk -v refpath="${rp}" '
BEGIN {
	OFS = "\t"
}
/^W/ { 
    if (refpath ~ $2) {
		print $5, $6, $7
	} else {
		print "Error: Expected path ID to match " refpath " but found " $2 " in file " FILENAME "\n" $0 > "/dev/stderr"
    	exit 1
	}
}
' tmp_${PREFIX}/chunk*.gfa > ${PREFIX}.walks.txt
rm -rf tmp_${PREFIX}


# Find the location of each fragment in the reference path (now assuming each barcode is unique, but if not only extra work, not error)
gawk -v prefix="${PREFIX}" '
BEGIN {
	OFS = "\t"
	output_file = prefix".fragment_positions.txt"
	errors_file = prefix".nodes_not_in_ref.txt"
}
# Function to get the ref position (same for F and R reads)
function ref_pos(read_aln) {
	delete aln_parts
	delete nodes
	split(read_aln, aln_parts, /_/)

	split(aln_parts[2], baseno, "-")
	dist_into_node = baseno[1]

	# Find whether the first node is in the reference path
	# If not, find the first node that is but note this had to be done. This also means using the first (or last)
	# base of the node rather than going by the position from the alignment
	split(aln_parts[1], nodes_only, /[<>]/)	# first element will be blank
	startnode = ""
	notfirst = 0
	if(pos[nodes_only[2]]) {
		startindex = 1
		startnode = nodes_only[2]
	} else {
		notfirst = 1
		dist_into_node = 0
		for (i = 3; i <= length(nodes_only); i++) {
			if (pos[nodes_only[i]]) {
				startindex = i-1
				startnode = nodes_only[i]
				break
			}
		}
	}
	if (!startnode) {
		return ""
	}

	# To get directions as well, replace each node with itself + colon separator
	# Then find the first node, check if direction is the same, and work out where start is
	# Since looking up by node ID, can check if node+direction is identical to field 3 from walks file
	# (if yes then same direction, if not then opposite direction)
    gsub(/[<>][0-9]+/, "&:", aln_parts[1])
    split(aln_parts[1], nodes, ":")
	split(pos[startnode], node_details)
	refpos = ""
	if (node_details[3] == nodes[startindex]) {
		# same direction so use start plus aln length
		refpos = node_details[1] + dist_into_node + 1	# +1 to make 1-based
	} else {
		# opposite direction so use end minus aln length
		refpos = node_details[2] - dist_into_node	# end is already 1-based (BED-like)
	}
	refpos = refpos "_" notfirst
	return refpos
}

# Load the walks
FNR == NR {
	nodeid = $3
	gsub(/^[<>]/, "", nodeid)
	pos[nodeid] = $0
	next
}
# Print the header for the new file
FNR == 1 && NR > FNR {
	print "Barcode", "Start", "End", "RefStart", "RefEnd", "RefStrand" > output_file
}
# Look up the location of the barcodes
{
	# Do F read (col 2)
	fpos = ref_pos($2)
	if (length(fpos) == 0) {
		print $1, "M" > errors_file
		next
	}

	# Do R read (col 3)
	rpos = ref_pos($3)
	if (length(rpos) == 0) {
		print $1, "M" > errors_file
		next
	}

	split(fpos, fres, "_" )
	split(rpos, rres, "_" )
	if (fres[2] || rres[2]) {
		print $1, "S" > errors_file
	}

	# Make start always less than end and give strand
	strand = "+"
	if (fres[1] > rres[1]) {
		temp = fres[1]
		fres[1] = rres[1]
		rres[1] = temp
		strand = "-"
	}
	print $1, $2, $3, fres[1], rres[1], strand > output_file

}

' ${PREFIX}.walks.txt ${PREFIX}.bc_to_check.txt

