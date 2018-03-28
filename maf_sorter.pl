#!/usr/bin/perl

# A simple script to sort MAF alignments by score.
# Score is simple: +1 to matches, -1 to mismatches and gaps.
# All alignments are put in a hash, and we make a score hash
# Sort hash by score, and then just grab alignments and print them properly

$align_counter = 1;

my %seq;
my %score;

open (READ, $ARGV[0]);
while (<READ>) {
    if (/^a\sscore/i) {
	# Found an alignment!
	my $line1 = <READ>;
	my $line2 = <READ>;
	my $score = &score($line1, $line2);

	my $alignment = $_;
	$alignment .= $line1;
	$alignment .= $line2;

	$seq{$align_counter} = $alignment;
	$score{$align_counter} = $score;
	$align_counter++;
    }
}
close READ;

for my $n (sort {$score{$b} <=> $score{$a}} keys %seq) {
    print $seq{$n}, "\n";
}

sub score {
    my ($s1, $s2) = @_;

    chomp ($s1);
    chomp ($s2);
    my @x = split (/\s+/, $s1);
    $s1 = $x[6];
    @x = split (/\s+/, $s2);
    $s2 = $x[6];

    my $total;
    $match = $s1 ^ $s2;
    while ($mask =~ /[^\0]/g) {
	$total++;
    }

    my $score = length($s1) - 4*$total;

    return $score;
}



