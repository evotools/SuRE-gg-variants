#!/usr/bin/env perl

use strict;
use warnings;

# impute_from_countstable.pl

# Go throught the fragments, check if any homozygous variants can be imputed (all) or if any heterozygous variants can be imputed (if haplotype already called from other variants)
# If a node is part of a variant, check the full path contains the full variant path

# To do imputation, try to find a match for the reads in the haplotype path. Don't do imputation if: reads don't match anywhere, if reads match more than one walk from a haplotype, if fragment is A but reads match only one haplotype, or if fragment is assigned as haplotype 1/2 but reads match haplotype 2/1 (respectively). However in these cases don't make any change to the fragment row, keep as it was before (including haplotype, alleles etc)

# Variants called in this way will have either an I in the caveats field (if other variants on the fragment, i.e. a kind of check the fragment is ok) or an II (no variants within the sequenced part).

# For multiallelic variants, search the fragment path for the unique paths as before (in which case they are marked I). If they aren't found then add the fragment with the allele inferred from the haplotype but mark A to indicate no sequence unique to this allele.

# As for vars_from_countstable.pl, use outputs from parse_var_path.pl and parse_var_path.multiallelic.pl

my $paths = shift;
my $fragments = shift;
my $vars = shift;
my $multi = shift;
my $colno = shift;	# 1-based column to find the genotype for this sample
my $out = shift;

my $coli = $colno - 1;

print "Checking for genotypes in var files column $colno\n\n";

open (OUT, ">$out") or die;
open (OUTS, ">fragstats.$out") or die;	# for recording lengths of fragments imputed from (first occurrence only)
print OUTS "Fpath\tRpath\tMaxNodes\tNvars\n";

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

# Load in variants in two hashes
my %vars_het;
my %vars_hom;
open (V, $vars) or die;
VAR:
while (<V>)	{
	next if ($_ =~ /RefChr/);
	my @ar = split (/\s+/, $_);
	my @paths = split (',', $ar[5]);
	my @gt = split (/\|/, $ar[$coli]);
	unless (defined $gt[0])	{
		die "Problem getting genotype from col index $coli:\n$_\n";
	}
	
	if ($ar[3] eq 'indel')	{
		my @limits = split (/\:/, $ar[4]);
		if ($gt[0] eq $gt[1])	{
			for (my $node = $limits[0]; $node <= $limits[1]; $node ++)	{
				$vars_hom{$node}{$ar[1]} = [($paths[$gt[0]], $gt[0], $ar[2], 'indel')];
			}
		}
		else	{
			for (my $node = $limits[0]; $node <= $limits[1]; $node ++)	{
				$vars_het{$node}{$ar[1]}{"1"} = [($paths[$gt[0]], $gt[0], $ar[2], 'indel')];
				$vars_het{$node}{$ar[1]}{"2"} = [($paths[$gt[1]], $gt[1], $ar[2], 'indel')];
			}
		}
	}
	elsif ($ar[3] eq 'snp' || $ar[3] eq 'mnp')	{
		my @snpnodes;
		for (my $i = 0; $i < @paths; $i ++)	{
			my @nodes = split (/(?=[<>])/, $paths[$i]);
			unless (@nodes == 3)	{
				warn "Too long MNP, skipping: $_\n";
				next VAR;
			}
			$nodes[1] =~ s/[<>]//;
			$snpnodes[$i] = $nodes[1];
		}
		if ($gt[0] eq $gt[1])	{
			$vars_hom{$snpnodes[$gt[0]]}{$ar[1]} = [($paths[$gt[0]], $gt[0], $ar[2], $ar[3])];
		}
		else	{
			$vars_het{$snpnodes[$gt[0]]}{$ar[1]}{"1"} = [($paths[$gt[0]], $gt[0], $ar[2], $ar[3])];
			$vars_het{$snpnodes[$gt[1]]}{$ar[1]}{"2"} = [($paths[$gt[1]], $gt[1], $ar[2], $ar[3])];
		}
	}
	elsif ($ar[3] eq 'uneqsub')	{
		my @nodes = split (/[\,<>]/, $ar[3]);
		if ($gt[0] eq $gt[1])	{
			foreach my $node (@nodes)	{
				next unless ($node);
				$vars_hom{$node}{$ar[1]} = [($paths[$gt[0]], $gt[0], $ar[2], $ar[3])];
			}
		}
		else	{
			foreach my $node (@nodes)	{
				next unless ($node);
				$vars_het{$node}{$ar[1]}{"1"} = [($paths[$gt[0]], $gt[0], $ar[2], $ar[3])];
				$vars_het{$node}{$ar[1]}{"2"} = [($paths[$gt[1]], $gt[1], $ar[2], $ar[3])];
			}
		}
	}
	else	{
		die "Problem: unknown variant type $ar[3] from $_\n";
	}

}

# Add non-biallelic variants to similar hashes
my %multi_het;
my %multi_hom;
open (M, $multi) or die;
MULTI:
while (<M>)	{
	next if ($_ =~ /RefChr/);
	my @ar = split (/\s+/, $_);
	my @paths = split (',', $ar[5]);
	my @gt = split (/\|/, $ar[$coli]);
	unless (defined $gt[0])	{
		die "Problem getting genotype from col index $coli:\n$_\n";
	}
	
	my @limits = split (/\:/, $ar[4]);
	my @flanks = split (/\:/, $ar[3]);
	if ($gt[0] eq $gt[1])	{
		for (my $node = $limits[0]; $node <= $limits[1]; $node ++)	{
			$multi_hom{$node}{$ar[1]} = [($paths[$gt[0]], $gt[0], $ar[2], 'multi', $flanks[0], $flanks[-1])];
		}
	}
	else	{
		for (my $node = $limits[0]; $node <= $limits[1]; $node ++)	{
			$multi_het{$node}{$ar[1]}{"1"} = [($paths[$gt[0]], $gt[0], $ar[2], 'multi', $flanks[0], $flanks[-1])];
			$multi_het{$node}{$ar[1]}{"2"} = [($paths[$gt[1]], $gt[1], $ar[2], 'multi', $flanks[0], $flanks[-1])];
		}
	}
}


# Go through fragments, collect consecutive lines with the same read paths then process once for each
open (FRAGS, $fragments) or die "Could not open $fragments";
my @fragdetails;	# collect fragment-specific outputs for fragments with the same read nodes (var details and SuRE data)
my $groupF;
my $groupR;
my @groupvars;		# haplotype, variant IDs, alleles, positions, types, caveats
my $groupvars_string;
my $groupscount = 0;

open (FRAGS, $fragments) or die;
FRAGMENT:
while (<FRAGS>)	{
	if ($_ =~ /barcode/i)	{
		print OUT $_;
		next;
	}
	my @frag = split (/\s+/, $_);
	# Check if read paths match current
	my $thisF = $frag[2] =~ s/_.*//r;
	my $thisR = $frag[3] =~ s/_.*//r;
	
	my @fragres = (join ("\t", @frag[0..6]), join ("\t", splice(@frag, 13)));
	
	my $thisvars = join ("\t", @frag[7..12]);
	
	# Do first assignation
	unless ($groupF)	{
		$groupF = $thisF;
		$groupR = $thisR;
		$groupvars_string = $thisvars;
		@groupvars = @frag[7..12];
		push (@fragdetails, \@fragres);
		next FRAGMENT;
	}
	
	if (($groupF eq $thisF) && ($groupR eq $thisR))	{
		unless ($thisvars eq $groupvars_string)	{
			die "Variant details in input should match but do not. Variants for other reads with this path: $groupvars_string\n$_";
		}
		push (@fragdetails, \@fragres);
		next FRAGMENT;
	}
	
	# If not, process last batch and assign group variables for next time
	my $varcols = process_group($groupF, $groupR, \@groupvars);
	for (my $i = 0; $i < @fragdetails; $i ++)	{
		print OUT "${$fragdetails[$i]}[0]\t$varcols\t${$fragdetails[$i]}[1]\n";
	}
	$groupscount ++;

	$groupF = $thisF;
	$groupR = $thisR;
	$groupvars_string = $thisvars;
	@groupvars = @frag[7..12];
	@fragdetails = (\@fragres);
	
}

# process last group
my $varcols = process_group($groupF, $groupR, \@groupvars);
for (my $i = 0; $i < @fragdetails; $i ++)	{
	print OUT "${$fragdetails[$i]}[0]\t$varcols\t${$fragdetails[$i]}[1]\n";
}
$groupscount ++;

warn "Number of groups: $groupscount\n";

# Process one pair of F and R read paths (nodes only)
sub process_group	{
	my ($F, $R) = @_;
	my $Frc = revcomp_path ($F);
	my $Rrc = revcomp_path ($R);
	
	# If reads meet in the middle then no missing sequence so no point trying to impute
	if ($F eq $Rrc)	{
		return ($groupvars_string);
	}
	my $meet = path_overlap($F, $Rrc);
	if ($meet)	{
		return ($groupvars_string);
	}
	
	# Don't impute one of the variants disagrees with the genotype (i.e. has caveat 'E')
	if ($groupvars[5] =~ /E/)	{
		return ($groupvars_string);
	}
	
	# If reads don't meet in the middle and don't conflict with the genotype, try to do imputation
	# First get possible full paths for the read and check ok for imputation
	my %matches;
	my %matches_rc;
	my $maxnodes = 0;
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
			return ($groupvars_string);
		}
		if ($hapmatches[0])	{
			$matches{$hap} = $hapmatches[0];
			$matches_rc{$hap} = revcomp_path($hapmatches[0]);
			my $count = $hapmatches[0] =~ tr/<>//;
			if ($count > $maxnodes)	{
				$maxnodes = $count;
			}
		}
	}
	
	# Check haplotypes matched are consistent with assignment from the reads, and get the nodes
	# For A get all nodes found in both haplotypes, and later check the variant path is found in both
	my $assigned = $groupvars[0];
	my @fragnodes;
	if ($assigned eq 'A')	{
		unless (scalar keys %matches == 2)	{
			return ($groupvars_string);
		}
		my %hap1;
		foreach my $node (split (/(?=[<>])/, $matches{"1"}))	{
			$node =~ s/[<>]//;
			$hap1{$node} = 1;
		}
		
		if ($matches{"1"} eq $matches{"2"})	{
			@fragnodes = keys %hap1;
		}
		else	{
			foreach my $node2 (split (/(?=[<>])/, $matches{"2"}))	{
				$node2 =~ s/[<>]//;
				if ($hap1{$node2})	{
					push (@fragnodes, $node2);
				}
			}
		}
	}
	elsif (($assigned eq "1") || ($assigned eq "2"))	{
		my $notassigned = 3 - $assigned;
		unless ($matches{$assigned} && (!$matches{$notassigned}))	{
			return ($groupvars_string);
		}
		foreach my $node (split (/(?=[<>])/, $matches{$assigned}))	{
			$node =~ s/[<>]//;
			push (@fragnodes, $node);
		}
	}
	else	{
		die "Problem with haplotype: @groupvars\n";
	}
	
	# Get variants already found
	my %got;
	foreach my $gotvar (split (/,/, $groupvars[1]))	{
		$got{$gotvar} = 1;
	}
	
	# Find expected variant path and check that it is present in all available paths
	# Since the walk may be in either direction relative to the ref path, need to check the variant path in both orientations
	# For homozygous variants can impute whether or not fragment allele has been assigned
	# For heterozygous variants, only check if haplotype is assigned
	# For multiallelic variants, possible for node to be part of the variant but have no unique paths. Want to know about this unless the fragment just contains one of the flanking nodes. If node is not a flanking node, add var to the fragment but with an A caveat. If node is a flanking one, just move to the next node -- if it's not just this node then var will be added anyway.
	# Variant is added to 'checked' for all biallelic. For multi it's mostly added but not if 1) there is no match and 2) the node is a flanking one, because in that case I want to check whether there are any non-flanking nodes as well
	my %new;
	my %checked;
	foreach my $node (@fragnodes)	{
		if ($vars_hom{$node})	{
			foreach my $varid (keys %{$vars_hom{$node}})	{
				next if ($got{$varid} || $checked{$varid});
				$checked{$varid} = 1;
				my $varpath = ${$vars_hom{$node}{$varid}}[0];
				my $matchcount = 0;
				foreach my $matchedhap (keys %matches)	{
					if (($matches{$matchedhap} =~ /$varpath(?=[<>]|\z)/) || ($matches_rc{$matchedhap} =~ /$varpath(?=[<>]|\z)/))	{
						$matchcount ++;
					}
				}
				if ($matchcount)	{
					if ($matchcount == (scalar keys %matches))	{
						$new{$varid} = [@{$vars_hom{$node}{$varid}}[1..3]];
					}
					else	{
						warn "Problem with matching path for homozygous variant (biallelic)\n$_\t$varid\t$varpath\n"; # since either both haplotypes (if there) should match or none
					}
				}
			}
		}
		if ($multi_hom{$node})	{
			foreach my $varid (keys %{$multi_hom{$node}})	{
				next if ($got{$varid} || $checked{$varid});
				my $varpath = ${$multi_hom{$node}{$varid}}[0];
				my $matchcount = 0;
				foreach my $matchedhap (keys %matches)	{
					if (($matches{$matchedhap} =~ /$varpath(?=[<>]|\z)/) || ($matches_rc{$matchedhap} =~ /$varpath(?=[<>]|\z)/))	{
						$matchcount ++;
					}
				}
				if ($matchcount == 0)	{
					unless ($node eq ${$multi_hom{$node}{$varid}}[4] || $node eq ${$multi_hom{$node}{$varid}}[5])	{
						$new{$varid} = [(@{$multi_hom{$node}{$varid}}[1..3], "A")];
						$checked{$varid} = 1;
					}
				}
				elsif ($matchcount == (scalar keys %matches))	{	# i.e. all haplotypes matched by reads also match the allele
					$new{$varid} = [@{$multi_hom{$node}{$varid}}[1..3]];
					$checked{$varid} = 1;
				}
				else	{
					warn "Problem with matching path for homozygous variant (multi)\n$_$varid\t$varpath\n";	# since either both haplotypes (if there) should match or none
					$checked{$varid} = 1;
				}
			}
		}
		if (($assigned eq "1") || ($assigned eq "2"))	{
			if ($vars_het{$node})	{
				foreach my $varid (keys %{$vars_het{$node}})	{
					next if ($got{$varid} || $checked{$varid});
					$checked{$varid} = 1;
					my $varpath = ${$vars_het{$node}{$varid}{$assigned}}[0];
					if (($matches{$assigned} =~ /$varpath(?=[<>]|\z)/) || ($matches_rc{$assigned} =~ /$varpath(?=[<>]|\z)/))	{
						$new{$varid} = [@{$vars_het{$node}{$varid}{$assigned}}[1..3]];
					}
				}
			}
			if ($multi_het{$node})	{
				foreach my $varid (keys %{$multi_het{$node}})	{
					next if ($got{$varid} || $checked{$varid});
					my $varpath = ${$multi_het{$node}{$varid}{$assigned}}[0];
					if (($matches{$assigned} =~ /$varpath(?=[<>]|\z)/) || ($matches_rc{$assigned} =~ /$varpath(?=[<>]|\z)/))	{
						$new{$varid} = [@{$multi_het{$node}{$varid}{$assigned}}[1..3]];
						$checked{$varid} = 1;
					}
					else	{
						unless ($node eq ${$multi_het{$node}{$varid}{$assigned}}[4] || $node eq ${$multi_het{$node}{$varid}{$assigned}}[5])	{
							$new{$varid} = [(@{$multi_het{$node}{$varid}{$assigned}}[1..3], "A")];
							$checked{$varid} = 1;
						}
					}
				}
			}
		}
	}

	# If any new variants were found, add to the variant details lines
	# column each for IDs, alleles, positions, types, caveats
	unless (keys %new)	{
		return ($groupvars_string);
	}
	my @newdetails;
	my $newcav;
	my $sep = ",";
	if ($groupvars[1] eq 'NA')	{
		@newdetails = ("") x 5;
		$sep = "";
		$newcav = "II";
	}
	else	{
		@newdetails = @groupvars[1..5];
		$newcav = "I";
	}
	
	foreach my $newvar (keys %new)	{
		$newdetails[0] .= $sep . $newvar;
		$newdetails[1] .= $sep . ${$new{$newvar}}[0];
		$newdetails[2] .= $sep . ${$new{$newvar}}[1];
		$newdetails[3] .= $sep . ${$new{$newvar}}[2];
		if (${$new{$newvar}}[3])	{
			$newdetails[4] .= $sep . "A";
		}
		else	{
			$newdetails[4] .= $sep . $newcav;
		}
		$sep = ",";
	}
	my $newres = "$groupvars[0]\t" . join ("\t", @newdetails);
	
	print OUTS "$groupF\t$groupR\t$maxnodes\t", (scalar (keys %new) + scalar (keys %got)), "\n";
	
	return($newres);
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