#!/usr/bin/perl

use strict;
use warnings;

my $usage = "$0 reference_genome list_of_maf_alignments\n";
our $word = 79;

my $output_chunks = "";
# ^ ^ ^
# We'll store the extra alignments here, in FASTA
# format, and print it at the end.
# Header will be: Sequence_Reference_position_CIGAR. Should be unique!


# List of MAF alignments should have:
#
# Whatever.whatever_something.maf Genome_ID

unless ($#ARGV == 1) {
    die $usage;
}

my %seq;
my $name;

open (READ, $ARGV[0]);
while (<READ>) {
    chomp;
    if (s/>//) {
	my @x = split;
	$name = $x[0];
	print STDERR $name;
    }
    else {
	$seq{$name} .= $_;
    }
}
close READ;
    
my @files;
my @genomes;
open (READ, $ARGV[1]);
while (<READ>) {
    chomp;
    my @x = split;
    push(@files, $x[0]);
    push(@genomes, $x[1]);
}
close READ;

# First pass, we'll be creating the table of variants
# and selecting which ones we will be adding to the
# reference genome.
#
# At the end of the day the proof of concept doesn't have
# to be perfect, if the hypothesis is sound then there
# will be improvements even if it misses stuff

for my $n (0..$#files) {

    # We're going to create a couple of hashes to track the intervals. Why?
    # Because that way there's a lot less to go through!

    my %ref_intervals;
    my %query_intervals;

    # For each sequence we'll create an array, and populate it initially with
    # the empty element (NONE 0 0);
    
    open (READ, $files[$n]);
    my %loopref = %seq;
    LOOP: while (<READ>) {
	if (/^a\sscore/i) {
	    $_ = <READ>;
	    chomp;
	    my @x = split;
	    my $ref = $x[$#x];

	    $_ = <READ>;
	    my @y = split;
	    my $query = $y[$#y];

	    if ($x[3] < 500) {
		next LOOP;
	    }

	    unless (exists $ref_intervals{$x[1]}) {
		my @anyny = ("NONE 0 0");
		$ref_intervals{$x[1]} = \@anyny;
	    }

	    # Intervals are too slow for the reference! Let's try another
	    # approach;

	    if (substr($loopref{$x[1]},$x[2],$x[3]) =~ /X+/) {
		next LOOP;
	    }
	    else {
		my $replacement = substr($loopref{$x[1]},$x[2],$x[3]);
		$replacement =~ s/-//g;
		$replacement =~ s/./X/g;
		substr($loopref{$x[1]},$x[2],length($replacement),$replacement);
	    }
	    
	    unless (exists $query_intervals{$y[1]}) {
		my @anyny = ("NONE 0 0");
		$query_intervals{$x[1]} = \@anyny;
	    }

	    $ref = uc($ref);
	    $query = uc($query);
	    
	    # Now we go through the intervals and check if any of the matches
	    # have already been used in one way or another! If they have,
	    # we discard this thing and continue

	    my $rarray = $ref_intervals{$x[1]};
	    my $qarray = $query_intervals{$y[1]};

	    # We're using the substring approach for the reference for the
	    # sake of speed; And we're not doing substrings for now.
	    
#	    foreach my $i (@$rarray) {
#		my @z = split(/\s/, $i);
#		    if ($z[1] eq $x[0]) {
#			if ($z[1] < $x[3] && $x[3] < $z[2]) {
#			    next LOOP;
#			}
#		}
#	    }
	    foreach my $i (@$qarray) {
		my @z = split(/\s/, $i);
		if ($z[0] eq $x[1]) {
		    if ($z[0] < $y[3] && $y[3] < $z[1]) {
			next LOOP;
		    }
		}
	    }


	    # So, not already explored! Let's go through the variants. 
	    print STDERR ".";

	    my $positions = length($ref);
	    my @SNPs;
	    my $finder = $ref ^ $query; # We're gonna XOR the sucker
	    while ($finder =~ m/[^\0]/g) {
		push (@SNPs, $-[0]);
		#print STDERR $-[0], "\n";
	    } # Easy! We have the positions of all the differences.
	    
	    # Next, we're going to get all stretches that have less than 80
	    # bp between them

	    my $start_slice = 0;
	    my $end_slice = 0;

	    for my $i (1..$#SNPs) {
		if ($SNPs[$i] - $SNPs[$end_slice] > $word || $i == $#SNPs) {

		    #print STDERR "i = $SNPs[$i]\nend = $SNPs[$end_slice]\nstr = $SNPs[$start_slice]\n";
		    
		    # Last slice is self-contained. We should export it!
		    # That means extract it, make a CIGAR, and call it a day
		    # Additionally, one last check if we're on the last SNP
		    # If that's the case, we just add an extra operation to
		    # include that last SNP. This avoids having to redo the
		    # whole damn thing when we leave the loop.

		    if ($i == $#SNPs) {
			$end_slice = $i;
		    }

		    my $sst;
		    my $lgth;

		    if ($SNPs[$start_slice] < $word) {
			$sst = 0;
		    }
		    else {
			$sst = $SNPs[$start_slice] - $word;
		    }

		    $lgth = $SNPs[$end_slice] - $sst + $word + 1;
		    
		    my $chunk = substr(
			$query,
			$sst,
			$lgth
			);
		    my $ref_chunk = substr(
			$ref,
			$sst,
			$lgth
			);
		    # Let's make the CIGAR;

		    my $cigar = &make_cigar ($ref_chunk, $chunk);

		    # Now we create the FASTA header

		    my $fasta_header = ">$y[1]_$x[1]_";
		    $fasta_header .= $x[2] + $SNPs[$start_slice] - $word;
		    $fasta_header .= "_$cigar\n";

		    $output_chunks .= $fasta_header;
		    $chunk =~ s/-//g;
		    $output_chunks .= "$chunk\n";

		    #Finally, add the intervals to the blacklist!

		    # Reference first
		    my $interval = "$x[1] ";
		    $interval .= $x[2] + $SNPs[$start_slice] - $word;
		    $interval .= " ";
		    $interval .= $x[2] + $SNPs[$end_slice] + $word;
		    push (@$rarray, $interval);

		    # Now the query, to avoid repeats/multimaps

		    $interval = "$y[1] ";
		    $interval .= $y[2] + $SNPs[$start_slice] - $word;
		    $interval .= " ";
		    $interval .= $y[2] + $SNPs[$end_slice] + $word;
		    push (@$qarray, $interval);

		    # And finally, we reset the intervals!
		    $start_slice = $i;
		    $end_slice = $i;
		}
		else {
		    $end_slice = $i;
		}
	    }
	}
    }
}

print STDERR "\n";
print $output_chunks;

sub make_cigar {
    my ($s1, $s2) = @_;
    my $cigarette = "";

    my $CIGAR = "";
    
    for my $i (0..length($s1)-1) {
	if (substr($s1,$i,1) ne "-" && substr($s2,$i,1) ne "-") {
	    $cigarette .= "M";
	}
	elsif (substr($s1,$i,1) eq "-" && substr($s2,$i,1) ne "-") {
	    $cigarette .= "I";
	}
	elsif (substr($s1,$i,1) ne "-" && substr($s2,$i,1) eq "-") {
	    $cigarette .= "D";
	}
    }

    while (length($cigarette) > 0) {
	my $symbol = substr($cigarette,0,1);
	$cigarette =~ s/^($symbol+)//g;
	$CIGAR .= length($1);
	$CIGAR .= $symbol;
    }

    return $CIGAR;
}
