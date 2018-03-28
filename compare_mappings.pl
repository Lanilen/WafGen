#!/usr/bin/perl

use strict;
use warnings;

my $usage = "$0 Reference sam1 sam2 [sam3 ...]

This program will read two or more SAM files and produce a series of
statistics about the differences between the mappings. It expects 'sam1' to be
mapped against the 'Reference' genome, and subsequent SAMs to be mapped
against additional genomes incorporated with the 'create_index.pl' script.

It will look for differences in mapping score for a given read, differences
in mapping position, and whether it has been mapped to 'sam2'+ only or has
been mapped against the original genome too.

Note that it will only consider mappings against sequences that are a single
mapping away from the original 'Reference' (so that they follow the convention
CONTIG_REF_POSITION_CIGAR). This behaviour is not expected to change, so
the way reference genomes are produced should reflect this.\n\n";

if ($#ARGV < 2) {
    die $usage;
}

open (READ, $ARGV[0]) || die;

my %reference;
my $name;
while (<READ>) {
    if (s/>//) {
	my @x = split;
	$reference{$x[1]} = 1;
    }
}
close READ;

my %original_score;
my %best_score;

my %original_mapping;
my %best_mapping;

my $original_gaps;
my $best_gaps;

my @new_mappings;

open (READ, $ARGV[1]) || die;
while (<READ>) {
    unless (/^\@/) {
	my @x = split;
	foreach my $n (@x) {
	    if ($n =~ /^AS\:/) {
		$n =~ /(\d+)$/;
		$original_score{$x[0]} = $1;
		$best_score{$x[0]} = $1;

		$original_mapping{$x[0]} = "$x[2] $x[3] $x[5]";
		$best_mapping{$x[0]} = "$x[2] $x[3] $x[5]";
	    }
	}
	while ($x[5] =~ m/(\d+)[ID]/g) {
	    $original_gaps += $1;
	}
    }
}
close READ;
$best_gaps = $original_gaps;

for my $i (2..$#ARGV) {
    my $current_gaps;
    open (READ, $ARGV[$i]) || die;
    while (<READ>) {
	unless (/^\@/) {
	    my @x = split;
	    foreach my $n (@x) {
		if ($n =~ /^AS\:/) {
		    $n =~ /(\d+)$/;
		    if (exists $best_score{$x[0]}) {
			if ($best_score{$x[0]} < $1) {
			    $best_score{$x[0]} = $1;
			    $best_mapping{$x[0]} = "$x[2] $x[3] $x[5]";
			}
		    }
		    else {
			$best_score{$x[0]} = $1;
		    }
		}
	    }
	    while ($x[5] =~ m/(\d+)[ID]/g) {
		$current_gaps += $1;
	    }
	}
    }
    close READ;
    if ($current_gaps < $original_gaps) {
	$best_gaps = $current_gaps;
    }
}

my $moved = 0;
my $changed = 0;
my $added = 0;
my $same = 0;
my $discarded = 0;
my $kept_score = 0;
my $better_score = 0;

foreach my $i (keys %best_mapping) {
    if (exists $original_mapping{$i}) {
	# There's both original and new mapping, so let's compare scores and
	# locations;
	if ($original_score{$i} == $best_score{$i}) {
	    $kept_score++;
	}
	elsif ($original_score{$i} < $best_score{$i}) {
	    $better_score++;
	}
	
	if ($original_mapping{$i} eq $best_mapping{$i}) {
	    $same++;
	}
	else {
	    my @x = split (/\s/, $best_mapping{$i});
	    # For reference, @x has:
	    # SEQID POSITION CIGAR
	    if ($x[0] =~ /\d+_\d+/) {
		my @trans = split (/_/, $x[0]);
		# We now have the "translation" from the mapping in the
		# reference sequence. If there's more than 4 elements it
		# is a cross-translation so we're going to skip it
		# for now; only consider the perfect 1:1 translations
		
		# For reference the ID of the translation should be:
		# CONTIG ORIGINALREF POSITION CIGAR
		# Where POSITION is the position the contig mapped to
		# the original reference.
		    
		if ($#trans == 3) {
		    # First, we translate the CIGAR into a string.
		    my $mapping;
		    while ($trans[3] =~ m/(\d+)([A-Z])/g) {
			$mapping .= ($2 x $1);
		    }
		    # Now we get a substr with the desired length being position of mapping -1
		    my $mapstring = substr($mapping,0,$x[1]-1);
		    my $originalpos = $trans[2];
		    $originalpos += ($mapstring =~ s/[MD]//g);
		    # Now check if the originalpos is within 10 bases of the original mapping
		    
		    my @original_coords = split (/\s/, $original_mapping{$i});
		    # We arbitrarily say 10 bp from the original mapping is the "same" mapping, because the ends
		    # of the reads tend to wobble a bit due to gap uncertainty.
		    
		    if ($original_coords[0] eq $trans[1] && abs($original_coords[1] - $originalpos) < 10) {
			$moved++;
		    }
		    else {
			$changed++;
		    }
		    
		}
		else {
		    $discarded++;
		}
	    }
	}
    }
    else {
	$added++;
    }
}

print "Moved:\t$moved
Changed:\t$changed
Added:\t$added
Same:\t$same
Discarded:\t$discarded\n
Kept Score:\t $kept_score
Better Score:\t$better_score
\n";
