#!/usr/bin/env perl

use strict;
use warnings;

$| = 1;

# find_frag_path.pl

# For fragments that map to significant variants, find the full fragment path so it can be plotted
# Edited down from impute_from_countstable.pl (so more complicated than it needs to be in places)
# Input is table from signif_frags.R of:
# Barcode, Start, End, Haplotype, library
# Output adds an extra column containing the assumed full fragment path, with a possible second extra column if there are two possibilities
# Also convert from path to comma-separated list of nodes
# Need to separate fragments by library before running and then run once for each library (since each will
# have different haplotype walks)

my $paths = shift;
my $fragments = shift;
my $out = shift;

open (OUT, ">$out") or die;

# Load in the two paths
my %paths;
my $sample;
open (P, $paths) or die "Could not open $paths";
while (<P>)	{
	my @ar  = split(/\s+/, $_);
	if ($sample)	{
		unless ($sample eq $ar[1])	{
			die "should be only one sample in walks file: $sample $ar[1]\n"
		}
	}
	else	{
		$sample = $ar[1];
	}
	my $hap = $ar[2];
	if ($paths{$hap})	{
		push (@{$paths{$hap}}, $ar[6]);
	}
	else	{
		${$paths{$hap}}[0] = $ar[6];
	}
}
close (P);

# Go through fragments, collect consecutive lines with the same read paths then process once for each
open (FRAGS, $fragments) or die "Could not open $fragments";
my @fragdetails;	# collect barcodes for fragments with the same read nodes (var details and SuRE data)
my $groupF;
my $groupR;
my $groupa;			# haplotype assignment
my $groupscount = 0;

open (FRAGS, $fragments) or die;
FRAGMENT:
while (<FRAGS>)	{
	if ($_ =~ /barcode/i)	{
		next;
	}
	my @frag = split (/\s+/, $_);
	# Check if read paths match current
	my $thisF = $frag[1] =~ s/_.*//r;
	my $thisR = $frag[2] =~ s/_.*//r;
	my $thisa = $frag[3];
	
	my $fragres = (join ("\t", @frag));
	
	# Do first assignment
	unless ($groupF)	{
		$groupF = $thisF;
		$groupR = $thisR;
		$groupa = $thisa;
		push (@fragdetails, $fragres);
		next FRAGMENT;
	}
	
	if (($groupF eq $thisF) && ($groupR eq $thisR))	{
		unless ($thisa eq $groupa) 	{
			die "Variant details in input should match but do not. $thisa $groupa\n$groupF $groupR\n$_";
		}
		push (@fragdetails, $fragres);
		next FRAGMENT;
	}
	
	# If not, process last batch and assign group variables for next time
	my $frag_path = process_group($groupF, $groupR);
#	print "Got $frag_path\n$groupF,$groupR\n";
	for (my $i = 0; $i < @fragdetails; $i ++)	{
		print OUT "$fragdetails[$i]\t$frag_path\n";
	}
	$groupscount ++;

	$groupF = $thisF;
	$groupR = $thisR;
	$groupa = $thisa;
	@fragdetails = ($fragres);
	
}

# process last group
my $frag_path = process_group($groupF, $groupR);
for (my $i = 0; $i < @fragdetails; $i ++)	{
	print OUT "$fragdetails[$i]\t$frag_path\n";
}
$groupscount ++;

warn "Number of groups: $groupscount\n";

# Process one pair of F and R read paths (nodes only)
sub process_group	{
	my ($F, $R) = @_;
	my $Frc = revcomp_path ($F);
	my $Rrc = revcomp_path ($R);
	
	my $sequenced;
	if ($F eq $Rrc)	{
		my $list = path2list($F);
		return ($list);
	}
	else	{
		my $meet = path_overlap($F, $Rrc);
		if ($meet)	{
			my $list = path2list($meet);
			return ($list);
		}
	}
	
	
	# If reads don't meet in the middle, try to get the full path of the fragment
	my %matches;
	foreach my $hap (sort keys %paths)	{
		my @hapmatches;
		foreach my $happath (@{$paths{$hap}})	{
			if ($happath =~ /($F[<>].*?$Rrc)(?=[<>]|\z)/ || $happath =~ /($R[<>].*?$Frc)(?=[<>]|\z)/)	{
				push (@hapmatches, $1);
			}
			elsif ($happath =~ /$F$Rrc(?=[<>]|\z)/)	{
				push (@hapmatches, "$F$Rrc");
			}
			elsif ($happath =~ /$R$Frc(?=[<>]|\z)/)	{
				push (@hapmatches, "$R$Frc");
			}
		}
		if (@hapmatches > 1)	{
			warn "Multiple matches in hap $hap for $F $R\n";
			return ("NA");
		}
		if ($hapmatches[0])	{
			$matches{$hap} = $hapmatches[0];
		}
	}
	
	# Check haplotypes matched are consistent with assignment from the reads
	my $assigned = $groupa;
	my @fragnodes;
	if ($assigned eq 'A')	{
		unless (scalar keys %matches == 2)	{
			return("NotMatched")
		}
		if ($matches{"1"} eq $matches{"2"})	{
			my $list = path2list($matches{"1"});
			return ($list);
		}
		else	{
			my $rced2 = revcomp_path($matches{"2"});
			if ($matches{"1"} eq $rced2)	{
				my $list = path2list($matches{"1"});
				return ($list);
			}
			else	{
				my $list = consensus($matches{"1"}, $matches{"2"});
		#	warn "done\n";
				return($list);
			}
		}
	}
	elsif (($assigned eq "1") || ($assigned eq "2"))	{
		my $notassigned = 3 - $assigned;
		unless ($matches{$assigned} && (!$matches{$notassigned}))	{
			return("NotMatched")
		}
		my $list = path2list($matches{$assigned});
	}
	else	{
		die "Problem with haplotype: $assigned $F $R\n";
	}
}


sub path_overlap	{
	my ($path1, $path2) = @_;
	
	my @ar1 = split(/(?=[<>])/, $path1);
	my @ar2 = split(/(?=[<>])/, $path2);

	my $length1 = scalar(@ar1);
	my $length2 = scalar(@ar2);
	my $len = $length2 - 1;
	my $off = $length1 - $length2;
	my $end = $length1 - 1;

	if ($off < 0)	{
		$off = 0;
		$len = $end;
	}
	while ($len >= 0)	{
		my $str1 = join('', @ar1[$off..$end]);
		my $str2 = join('', @ar2[0..$len]);
	
		if ($str1 eq $str2)	{
			# get the joined sequence (need to check)
			my $joinedpaths = $path1;
			my $ind = $len + 1;
			while ($ind < $length2)	{
				$joinedpaths .= $ar2[$ind];
				$ind ++;
			}
			return ($joinedpaths);
		}
		else	{
			$len --;
			$off ++;
		}
	}
	return(0);
}

sub revcomp_path	{
	(my $orig) = @_;
	my @nodes = (split(/(?=[<>])/, $orig));
	for (my $i = 0; $i < @nodes; $i ++)	{
		$nodes[$i] =~ tr/<>/></;
	}
	my $new = join ('', reverse(@nodes));
	return ($new);
}

sub path2list {
	(my $pathstr) = @_;
#	print "$pathstr\n";
	my @ar = split (/[<>]/, $pathstr);
#	print "@ar\n";
	shift @ar;
	my $liststr = join(',', @ar);
#	print "$liststr\n";
#	exit;
	return ($liststr);
}

# For cases where a fragment matches multiple non-identical haplotypes, get a consesnsus sequence where they go through the same nodes, and ? where they differ
sub consensus	{
	my ($path1, $path2) = @_;
	
	my @nodes1 = split (/[<>]/, $path1);
	shift (@nodes1);
	my @nodes2 = split (/[<>]/, $path2);
	shift (@nodes2);
	
	unless ($nodes1[0] eq $nodes2[0])	{
		if ($nodes1[0] eq $nodes2[-1])	{
			@nodes2 = reverse(@nodes2);
		}
		else	{
			die "Problem with paths: $path1\n$path2\n";
		}
	}

	my $i = 0;
	my $j = 0;
	my @comp;
	
#	print "$path1\n$path2\n";

	COMP:
	while ($i < @nodes1 && $j < @nodes2) {
		if ($nodes1[$i] eq $nodes2[$j]) {
			push (@comp, $nodes1[$i]);
			$i++;
			$j++;
			next;
		}
#		unless ($i > 0 & $comp[-1] eq '?')	{
    		push (@comp, '?');
# 	   }
		my $k = 1;
		while ($nodes1[$i+$k])	{
			if ($nodes1[$i+$k] eq $nodes2[$j])	{
				push (@comp, $nodes1[$i+$k]);
				$i = $i+$k + 1;
				$j ++;
				$k ++;
				next COMP;
			}
			else	{
				$k ++;
			}
		}
		$k = 1;
		while ($nodes2[$j+$k]) {
	#		unless ($nodes1[$i] && $nodes2[$j+$k])	{
	#			die "Problem: $i $j $k\n@nodes1\n@nodes2\n";
	#		}
			if ($nodes1[$i] eq $nodes2[$j+$k]) {
				push (@comp, $nodes1[$i]);
				$j = $j+$k;
				next COMP;
			}
			else	{
				$k ++;
			}
		}
		$i ++;
		$j ++;
	}
	
	my $cons = join(',', @comp);
	return($cons);
}

