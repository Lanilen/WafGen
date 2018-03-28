#!/usr/bin/perl

my @fastq = qw(SRR1011205_1.fastq.gz  SRR1011287_1.fastq.gz  SRR1210129_1.fastq.gz);


foreach my $n (@fastq) {
#    system("/home/roger/Downloads/NextGenMap-0.5.5/bin/ngm-0.5.5/ngm-core -t 14 -g 0,1 --max-read-length 80 --no-unal -r Reference.fna -q fastq/$n --very-fast -e -i 0.85 > /media/roger/Hanna/Plasmodium/0.85/0.$n");
    for my $i (1..5) {
#	system("/home/roger/Downloads/NextGenMap-0.5.5/bin/ngm-0.5.5/ngm-core -t 14 -g 0,1 --max-read-length 80 --no-unal -r plus$i/Reference.fna -q fastq/$n --very-fast -e -i 0.85 > /media/roger/Hanna/Plasmodium/0.85/$i.$n");
	system("/home/roger/Downloads/NextGenMap-0.5.5/bin/ngm-0.5.5/ngm-core -t 14 -g 0,1 --max-read-length 80 --no-unal -r plus$i/query -q fastq/$n --very-fast > /media/roger/Hanna/Plasmodium/fullRefs/$i.$n.sam");
    }
}
