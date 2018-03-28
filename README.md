# WafGen

A simplified approach to create a pan-genome reference.

It requires the use of an external whole genome aligner and to provide the results in MAF format.

create_index.pl: The main script to create the pan-genome from a MAF file and some FASTA files.
                 It looks for differences between the reference and the new genome through the
                 MAF alignment, and extract divergent regions. Each region will have the differences
                 flanked by enough conserved sequence to encompass a whole Illumina read. This
                 parameter can easily be changed.
                 
maf_sorter.pl: A helper script that will take a MAF alignment and sort it by score, so that the
               longer, higher scoring pairs are considered first. Each region of the two genomes
               can only be considered ONCE for alignment, so the order in which they are considered
               is up to the user. By default, most tools will use input order to create the MAF
               output.
               
compare_mappings.pl: Given a reference and a series of SAM files, this script will give some statistics
                     on how many reads map with the same score, different, map in a different region,
                     etc.
