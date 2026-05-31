#!/usr/bin/env perl

use strict;
use warnings;

# parse_var_path.multiallelic.pl

# Take all alleles from a variant with >1 alt allele and find a path or a set of paths that are unique to each one
# Filter vcf first to remove variants with missing calls and alleles that don't actually appear in the sample genotypes, then to get variants that still have at least 3 alleles (including ref) actually used in the genotypes
# These should now not include any variant where all haplotypes have the same alt allele (since then there would only be one alt left) but anyway put in a command to ensure no allele appears more than total haps-1 times
# bcftools view -g ^miss -a <input.vcf.gz> | bcftools view --min-alleles 3 --max-ac <total haps minus 1>:major | bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AT\t[%GT\t]\n' > <output.txt>

# For each allele get a pipe-separarated string of what's unique to that allele (e.g. ">1|>3>4" means node >1 only appears in this allele, and while 3 and 4 might appear in others they only appear next to each other in this allele). Each will be searched for in the read path to see if the fragment comes from that allele. The pipe separation means the whole thing can be used directly in a regular expression to search for them all at once.

# Check to see if node range consecutive? Try and see if many are excluded

# Output
# Chr	VarID	pos	flanking	noderange	<list of how to find each allele>	genotypes
# In the 'how to find' column the alleles are separated by commas, and within each allele the features are pipe separated (if >1 is unique to that allele)

my $vars = shift;
my $maxpath = shift;	# to set a maximum number of nodes to check the unique path for; 0 if no max
my $samplenames = shift;	# comma-separated list of sample names to use for heading - same order as vcf
my $outname = shift;

# do something to keep a record of paths where I've not found unique features
open (OUTW, ">no_unique.txt") or die;
open (OUT, ">variants.$outname.txt") or die;

$samplenames =~ s/\,/\t/g;
print OUT "RefChr\tID\tRefPos\tFlanking\tRange\tAllelePaths\t$samplenames\n";

open (V, $vars) or die;
VARIANT:
while (<V>)	{
	chomp;
	my @info = split(/\s+/, $_);
	my @paths = split(/\,/, $info[5]);
	my @altseqs = split(/\,/, $info[4]);
	my @seqs = ($info[3], @altseqs);
	my $genotypes = join("\t", splice(@info, 6));
	
	# Split up all and get/check flanking nodes
	my @pathsall;
	my $flank1;
	my $flank2;
	my $deletion; # if one of the paths is just the flanking nodes this gives the index of that allele
	# If only two nodes, check that this is the only path for which that's true. (Should check if they are the flanks but this will be done through processing all the paths)
	my @unique = ("") x scalar(@paths);	# string of what's unique to each allele
	my @allnodes;	# can have repeats
	for (my $i = 0; $i < @paths; $i ++)	{
		my @nodes = split (/(?=[<>])/, $paths[$i]);
		push (@allnodes, @nodes);
		$pathsall[$i] = [@nodes];
		if ($i == 0)	{
			$flank1 = $nodes[0];
			$flank2 = $nodes[-1];
		}
		elsif (!(($nodes[0] eq $flank1) && ($nodes[-1] eq $flank2))) {
			die "Problem with flanking nodes: $_";
		}
		if (@nodes == 2)	{
			if (defined ($deletion))	{
				die "Problem: multiple 2-node paths\n$_";
			}
			$deletion = $i;
			$unique[$i] = $paths[$i];
		}
	}
	
	my ($is_sequence, $start, $end) = check_node_sequence(@allnodes);
	unless ($is_sequence)	{
		warn "Nodes out of sequence: $_\n";
		next VARIANT;
	}
	
	my @nounique;
	# for each path with >2 nodes, find nodes or sub-paths unique to it
	for (my $i = 0; $i < @paths; $i ++)	{
		next if (defined ($deletion) && $i == $deletion);
	#	print "doing allele $i\n";
		my %othernodes;
		for (my $j = 0; $j < @paths; $j ++)	{
			next if ($i == $j);
			foreach my $othernode (@{$pathsall[$j]})	{
				$othernodes{$othernode} = 1;
			}
		}
		# check whether any (non-flank) nodes in this path is unique to it
		my @oneunique;
		my %included_ind;
		# Check for the shortest unique paths in each allele
		# Start with 1-node paths, then search for 2-node paths among nodes that aren't already included, then 3- node paths, etc
		my $length = 1;
		my $max_path_length = scalar(@{$pathsall[$i]});
		if ($maxpath && ($maxpath < $max_path_length))	{
			$max_path_length = $maxpath;
		}
		while (($length <= $max_path_length) && ((scalar(@{$pathsall[$i]})-(scalar keys %included_ind) >= $length)))	{
			my $last = scalar @{$pathsall[$i]} - $length;
			for (my $j = 0; $j <= $last; $j ++)	{
				my $included = 0;
				my $check = $j;
				my $endnode = $j+$length-1;
				while ($check <= $endnode)	{
					if ($included_ind{$check})	{
						$included = 1;
						last;
					}
					$check ++;
				}
				next if ($included);
				my $subpath = join ('', @{$pathsall[$i]}[$j..$endnode]);
				my $found = 0;
				for (my $k = 0; $k < @paths; $k ++)	{
					next if ($i == $k);
					if ($paths[$k] =~ /$subpath/)	{
					#	print "found hit for $subpath in $paths[$k]\n";
						$found = 1;
					}
				}
				unless ($found)	{
					push (@oneunique, $subpath);
					my $inc = $j;
					while ($inc <= $endnode)	{
						$included_ind{$inc} = 1;
						$inc ++;
					}
				}
			}
			$length ++;
		}
		if (@oneunique)	{
			$unique[$i] = join ("|", @oneunique);
		}
		else	{
			print OUTW "Allele$i\t$_\n";
			next VARIANT;
		}
	}
	
	# output
	my $flankcol = "$flank1:$flank2";
	$flankcol =~ s/[<>]//g;
	print OUT "$info[1]\t$info[0]\t$info[2]\t$flankcol\t$start:$end\t", join(",", @unique), "\t$genotypes\n";
}

close (V);
#close (OUT);
close (OUTW);

# Subroutine to check whether the set of nodes forms a continuous sequence, ignoring repeats. Doesn't check if they are in order.
sub check_node_sequence	{
	my (@all) = @_;
	foreach my $one (@all)	{
		$one =~ s/[<>]//;
	}
	my @list = sort {$a <=> $b} @all;
	my ($start, $end) = ($list[0], $list[-1]);
	my $fail = 0;
	my @missing;
	for (my $i = 0; $i < @list-1; $i ++)	{
		next if (($list[$i] + 1 == $list[$i+1]) || ($list[$i] == $list[$i+1]));
		my $miss = $list[$i] + 1;
		while ($miss < $list[$i+1])	{
			push (@missing, $miss);
			$miss ++;
		}
		last;
	}
	if ($missing[0])	{
		return (0)
	}
	else	{
		return (1, $start, $end);
	}
}
