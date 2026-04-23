#!/usr/bin/env perl

use strict;
use warnings;

# vars_from_countstable.pl
# Check the read paths from the counts table to find any variants and (where possible) check which allele

# For biallelic variants use the output from parse_var_path.pl
# For multiallelic variants use the output from parse_var_path.multiallelic.pl
# To set which sample is being processed, give the 1-based column number in the biallelic variants file

# Output should be the same as the final var_counts_tables files but just has not had indels or imputed variants added yet (and also includes fragments with no variants)
# This is the same as the counts table except that the variant details are added between the variant details columns and the SURE counts columns (I changed this relative to the original counts tables so that everything except the SuRE results would always be in the same place even if the replicates changed)
# Chr     Barcode Start   End     First_node      Strand  MAPQs   Haplotype       VarID   VarAllele       VarPos  VarType VarCav  iPCR   SuRE_replicate1	SuRE_replicate2...

# VarCav is variant caveats -- reasons not to include this fragment if filtering. This could be 0 if fragment fine, E if the allele from the read is inconsistent with what's expected from the sample (e.g. allele 1 in homozygous ref), I if the allele is imputed (but don't have any of them yet)
# If different variants on the same allele give conflicting haplotypes then the fragment will not be output (which means that indels on this fragment won't be matched up with anything in the next stage, so they won't be output either). However record which variants are involved in the conflict and output if they appear in warnings more than a certain number of times in the whole data set (as this might imply a problem with the genotype rather than a sequencing error). N.B. The way this is coded I might get different counts from different runs because after the first conflict is found no more nodes are processeed on that fragment, and the nodes might be processed in a different order between runs. So only approximate but if there is an incorrect genotype then it will be counted every time it appears on a fragment, while other vars on the fragment might sometimes be counted and sometimes not.
# Also record when a fragment has nodes involved with a non-biallelic variant but no unique path, as want to count how many times this happens for each variant to see if it's a problem. NOT DONE?

my $vars = shift;
my $multi = shift;
my $colno = shift;	# 1-based column to find the genotype for this sample
my $out = shift;
my $coli = $colno - 1;
my $conflict_report = 5;	# report a variant if it appears in haplotype conflicts more than this many times

print "Checking for genotypes in vars and multiallelic files column $colno\n\n";

open (OUT, ">$out") or die;
open (OUTERR, ">warnings.$out") or die;
open (OUTAMBIG, ">ambiguous.$out") or die;

# Make two hashes of variants, SNP/MNPs and indels/uneqsubs (alleles are different lengths but neither allele is simply the two flanking nodes joined together). Difference is that for SNPs/MNPs one node distinguishes between the alleles (skipping MNPs where this is not the case, which is hardly any), while for the others there is more than one
# For SNPs/MNPs, variant is added to the node key if that node varies between alleles, and it then points to the full path for that allele and the haplotype, so don't need to check reads that don't have the important node, or check for both paths when you know it definitely doesn't have one of them
# For indels/subs, any node in the variant points to a hash of the unique paths for each allele (i.e. if this path is in the read you can say it's from this allele)
# If homozygous, allele present gives 'A' (for ambiguous), allele that doesn't appear gives 'E' (for error)
# If heterozygous, each allele points to its haplotype number
open (V, $vars) or die;
my %snps;	# and mnps
my %indels;	# and uneqsubs
VAR:
while (<V>)	{
	next if ($_ =~ /RefChr/);
	my @arv = split(/\s+/, $_);
	my @paths = split (',', $arv[5]);
	unless (@paths == 2)	{
		die "Too many alleles: $_\n";
	}
	
	my @gt = split (/\|/, $arv[$coli]);
	unless (defined $gt[0])	{
		die "Problem getting genotype from col index $coli:\n$_\n";
	}
	
	if ($arv[3] eq 'snp' || $arv[3] eq 'mnp')	{
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
			$snps{$snpnodes[$gt[0]]}{$arv[1]} = [($paths[$gt[0]], "A", $gt[0], $arv[2], $arv[3])];
			my $other = 1-$gt[0];
			$snps{$snpnodes[$other]}{$arv[1]} = [($paths[$other], "E", $other, $arv[2], $arv[3])];
		}
		else	{
			$snps{$snpnodes[$gt[0]]}{$arv[1]} = [($paths[$gt[0]], "1", $gt[0], $arv[2], $arv[3])];
			$snps{$snpnodes[$gt[1]]}{$arv[1]} = [($paths[$gt[1]], "2", $gt[1], $arv[2], $arv[3])];
		}
	}
	elsif ($arv[3] eq 'uneqsub')	{
		my @nodes = split (/[\,<>]/, $arv[5]);
		my %nodes;
		foreach my $arrnode (@nodes)	{
			next unless ($arrnode);
			$nodes{$arrnode} = 1;
		}
	
		# Put all nodes in either path into a hash with haplotype information
		# node-VarID-path-(haplotype, hap allele, var position)
		if ($gt[0] eq $gt[1])	{
			foreach my $node (keys %nodes)	{
				$indels{$node}{$arv[1]}{$paths[$gt[0]]} = [("A", $gt[0], $arv[2], $arv[3])];
				my $other = 1-$gt[0];
				$indels{$node}{$arv[1]}{$paths[$other]} = [("E", $other, $arv[2], $arv[3])];
			}
		}
		else	{
			foreach my $node (keys %nodes)	{
				$indels{$node}{$arv[1]}{$paths[$gt[0]]} = [("1", $gt[0], $arv[2], $arv[3])];
				$indels{$node}{$arv[1]}{$paths[$gt[1]]} = [("2", $gt[1], $arv[2], $arv[3])];
			}
		}
	}
	
	elsif ($arv[3] eq 'indel')	{
		my @limits = split (/\:/, $arv[4]);
	
		if ($gt[0] eq $gt[1])	{
			for (my $node = $limits[0]; $node <= $limits[1]; $node ++)	{
				$indels{$node}{$arv[1]}{$paths[$gt[0]]} = [("A", $gt[0], $arv[2], $arv[3])];
				my $other = 1-$gt[0];
				$indels{$node}{$arv[1]}{$paths[$other]} = [("E", $other, $arv[2], $arv[3])];
			}
		}
		else	{
			for (my $node = $limits[0]; $node <= $limits[1]; $node ++)	{
				$indels{$node}{$arv[1]}{$paths[$gt[0]]} = [("1", $gt[0], $arv[2], $arv[3])];
				$indels{$node}{$arv[1]}{$paths[$gt[1]]} = [("2", $gt[1], $arv[2], $arv[3])];
			}
		}
	}
	else	{
		die "Problem: unknown variant type $arv[3] from $_\n";
	}
		
}


my %multi_nodes;	# each node points to all the variants that involve it
my %multi_vars;
# multi_vars is a hash of variant IDs
# each has an array of arrays
# first is an array of the variant paths in allele order
# then an array of the haplotypes for those alleles (1, 2, A (homozygous for this allele in sample), E (this sample shouldn't have this allele)
# then an array of the two flanking nodes
# then the variant position
open (M, $multi) or die;
while (<M>)	{
	my @arm = split(/\s+/, $_);
	
	my @limits = split (/\:/, $arm[4]);
	for (my $nodeid = $limits[0]; $nodeid <= $limits[1]; $nodeid ++)	{
		$multi_nodes{$nodeid}{$arm[1]} = 1;
	}
		
	# get the variant details
	my @paths = split(/\,/, $arm[5]);
	my @gt = split (/\|/, $arm[$coli]);
	unless (defined $gt[0])	{
		die "Problem getting genotype from col index $coli:\n$_\n";
	}
	
	my @haps;
	if ($gt[0] == $gt[1])	{
		$haps[$gt[0]] = "A";
	}
	else	{
		$haps[$gt[0]] = "1";
		$haps[$gt[1]] = "2";
	}
	for (my $i = 0; $i < @paths; $i ++)	{
		unless (defined $haps[$i])	{
			$haps[$i] = "E";
		}
	}
	my @flanks = split(/\:/, $arm[3]);
	$multi_vars{$arm[1]} = [([@paths],[@haps],@flanks, $arm[2])];
}
close (M);

# For each node in the fragment, check if it is associated with a variant and get the details
# If reads don't overlap then do separately and check any shared variants have the same allele (for long variants could get two non-overlapping reads covering different bits of it)
#open (FRAGS, $fragments) or die "Could not open $fragments";
my %conflicts;
FRAGMENT:
while (<STDIN>)	{
	my @frag = split (/\s+/, $_);
	if ($_ =~ /Barcode/i)	{
		print OUT join ("\t", @frag[0..6]), "\tHaplotype\tVarID\tVarAllele\tVarPos\tVarType\tVarCav\t", join ("\t", splice(@frag, 7)), "\n";
		next;
	}
	my $F = $frag[2] =~ s/_.*//r;
	my $R = $frag[3] =~ s/_.*//r;
	my $Rrc = revcomp_path ($R);
	
	my $resline = $frag[0];
	
	my @fragpaths;
	# Check if paths are identical or overlapping
	if ($F eq $Rrc)	{
		$fragpaths[0] = $F;
	}
	else	{
		my $meet = path_overlap($F, $Rrc);
		if ($meet)	{
			$fragpaths[0] = $meet;
		}
		else	{
			@fragpaths = ($F, $R);
		}
	}
	
	# Get variants for each allele
	# Collect in %fragvars then can check as you go along if it's been found already and if the two versions match. Get array of haplotype, allele, position, type for each
	my %fragvars;
	my $hap;
	foreach my $fragpath (@fragpaths)	{
		my %pathvars;
		my $rc = revcomp_path($fragpath);
		my @nodes = split (/(?=[<>])/, $fragpath);
	#	print "nodes @nodes\n";
		my %multi_near;
		my %readvars;
		foreach my $node (@nodes)	{
			$node =~ s/[<>]//;
			if ($snps{$node})	{
				foreach my $varid (keys %{$snps{$node}})	{
					if ($readvars{$varid})	{
						warn "More than one SNP/MNP allele from single read: $_\n";
						next FRAGMENT;
					}
					my $varpath = ${$snps{$node}{$varid}}[0];
					if (($fragpath =~ /$varpath(?=[<>]|\z)/) || ($rc =~ /$varpath(?=[<>]|\z)/))	{
						$readvars{$varid} = [@{$snps{$node}{$varid}}[1..4]];
					}
					else	{
						$readvars{$varid} = "nomatch";
					}
				}
			}
			if ($indels{$node})	{
				foreach my $varid (keys %{$indels{$node}})	{
					next if ($readvars{$varid});
					foreach my $varpath (keys %{$indels{$node}{$varid}})	{
						if (($fragpath =~ /$varpath(?=[<>]|\z)/) || ($rc =~ /$varpath(?=[<>]|\z)/))	{
							$readvars{$varid} = [@{$indels{$node}{$varid}{$varpath}}];
						}
					}
					unless ($readvars{$varid})	{
						$readvars{$varid} = "nomatch";
					}
				}
			}
			# Find nodes that are involved with multiallelic variants (process after all found)
			if ($multi_nodes{$node})	{
				foreach my $varid (keys %{$multi_nodes{$node}})	{
					if ($multi_near{$varid})	{
						push (@{$multi_near{$varid}}, $node);
					}
					else	{
						${$multi_near{$varid}}[0] = $node;
					}
				}
			}
		}
		
		# If any multiallelic variants involve these nodes, check the read path for unique paths for the alleles
		# If none are found, and more than just a flanking node is in the read, record that this fragment overlaps the variant but its genotype can't be determined (or it has a mismatch)
		if (keys %multi_near)	{
			foreach my $varid (keys %multi_near)	{
				for (my $i = 0; $i <@{$multi_vars{$varid}[0]}; $i ++)	{
					my $uniquepath = ${$multi_vars{$varid}[0]}[$i];
					if (($fragpath =~ /$uniquepath(?=[<>]|\z)/) || ($rc =~ /$uniquepath(?=[<>]|\z)/))	{
						$readvars{$varid} = [(${$multi_vars{$varid}[1]}[$i], $i, ${$multi_vars{$varid}}[4], "multi")];
						last;
					}
				}
				unless ($readvars{$varid})	{
					if (@{$multi_near{$varid}} > 1) 	{
						# fragment is from this variant but can't genotype
						print OUTAMBIG "Ambiguous1:\t", join ("\t", @frag[0..3]), "\t$varid\n";
					}
					else	{
						my $tomatch = ${$multi_near{$varid}}[0];
						$tomatch =~ s/[<>]//;
						unless ($tomatch eq ${$multi_vars{$varid}}[2] || $tomatch eq ${$multi_vars{$varid}}[3])	{
						# only one node but it's not a flanking node, so again fragment is from this variant but can't genotype
							print OUTAMBIG "Ambiguous1:\t", join ("\t", @frag[0..3]), "\t$varid\n";
						}
					}
				}
			}
		}
		# Check consistent with any variants found on the fragment so far
		# If same variant found on multiple reads check the same haplotype and allele
		# For all variants check for conflicts over haplotype (can be ambiguous or E but if assigned must agree)
		foreach my $varid (keys %readvars)	{
			next if ($readvars{$varid} eq 'nomatch');
			if ($fragvars{$varid})	{
				if ((${$readvars{$varid}}[1] ne ${$fragvars{$varid}}[1]) || (${$readvars{$varid}}[0] ne ${$fragvars{$varid}}[0]))	{
					print OUTERR "Pair conflict:\t", join (",", @{$readvars{$varid}}), "\t", join(",", @{$fragvars{$varid}}), "\n";
					next FRAGMENT;
				}
			}
			else	{
				if ((${$readvars{$varid}}[0] eq "1") || (${$readvars{$varid}}[0] eq "2"))	{
					if (!(defined $hap))	{
						$hap = ${$readvars{$varid}}[0];
					}
					elsif ($hap ne ${$readvars{$varid}}[0])	{
						print OUTERR "Haplotype conflict:\t", join ("\t", @frag[0..3]), "\t", join (",", keys %fragvars), ",$varid\t$hap,${$readvars{$varid}}[0]\n";
						foreach my $varid (keys %fragvars, $varid)	{
							if ($conflicts{$varid})	{
								$conflicts{$varid} ++
							}
							else	{
								$conflicts{$varid} = 1;
							}
						}
						next FRAGMENT;
					}
				}
				$fragvars{$varid} = $readvars{$varid};
			}
		}
	}
	# Output entire fragment (adding VarID   VarAllele       VarPos  VarType VarCav to original row)
	my @extra_cols;
	if (keys %fragvars)	{
		@extra_cols = "" x 5;
		my @varids = sort keys %fragvars;
		my $sep = ",";
		for (my $i = 0; $i < @varids; $i ++) {
			if ($i+1 == @varids)	{
				$sep = "";
			}
			$extra_cols[0] .= $varids[$i] . $sep;
			$extra_cols[1] .= ${$fragvars{$varids[$i]}}[1] . $sep;
			$extra_cols[2] .= ${$fragvars{$varids[$i]}}[2] . $sep;
			$extra_cols[3] .= ${$fragvars{$varids[$i]}}[3] . $sep;
			if (${$fragvars{$varids[$i]}}[0] eq "E")	{
				$extra_cols[4] .= "E" . $sep;
			}
			else	{
				$extra_cols[4] .= "0" . $sep;
			}
		}
	}
	else	{
		@extra_cols = ("NA") x 5;
	}
	unless (defined $hap)	{
		$hap = "A";
	}
	
	print OUT join ("\t", @frag[0..6]), "\t$hap\t", join("\t", @extra_cols), "\t";
	my @counts = splice (@frag, 7);
	print OUT join("\t", @counts), "\n";
}
close (OUT);
close (OUTERR);
close (OUTAMBIG);

open (C, ">conflicts.$out");
foreach my $conf (keys %conflicts)	{
	if ($conflicts{$conf} > $conflict_report)	{
		print C "$conf\t$conflicts{$conf}\n";
	}
}
close (C);

# check if they overlap and give the merged path
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
