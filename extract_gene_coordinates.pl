#!/usr/bin/perl
use warnings;
use strict;
use List::MoreUtils qw( minmax );

# Extract gene coordinates from CDS/exon GTF file

my (@F, @last, @positions);
my ($start, $end);
while (<>) {
    @F = split /\t/;
    next if /^#/ or not @last;
    my ($gene_F) = $F[8] =~ m/gene_id "(.*?)"/;
    my ($gene_last) = $last[8] =~ m/gene_id "(.*?)"/;
    push(@positions, $last[3], $last[4]);
    if ( $gene_last ne $gene_F ) {
        ($start, $end) = minmax @positions;
        splice @last, 2, 1, 'gene';
        splice @last, 3, 1, $start;
        splice @last, 4, 1, $end;
        splice @last, 8, 1, $gene_last;
        print join("\t", @last), "\n";
        splice( @positions ); 
    }
} continue {
    @last = @F;
}

push(@positions, $last[3], $last[4]);
my ($gene_last) = $last[8] =~ m/gene_id\s"(.*?)"/;
($start, $end) = minmax @positions;
splice @last, 2, 1, 'gene';
splice @last, 3, 1, $start;
splice @last, 4, 1, $end;
splice @last, 8, 1, $gene_last;
print join("\t", @last), "\n";
        
