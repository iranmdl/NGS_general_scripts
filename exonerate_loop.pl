#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;

my $genes = shift;
my $chr_path = "/media/sequentia/NAS/projects/042-melon_puig/exonerate/3.6_chr_single_files/";

sub SingleFasta {
   my $file = shift;
   open ( FILE, $file ) or die "Could not open file $file: $!";
   # change the Input Record Separator to ">" which allows to take in a full, multiline
   # FASTA record at once.
   $/ = ">";
   my $junk = <FILE>; # Discard the ">" at the begining of the file
   # Now read through your input file one sequence record at a time. Each input record will be a
   # multiline FASTA entry.
   while ( my $record = <FILE> ) {
      chomp $record; 
      my ($defLine, @seqLines) = split /\n/, $record;
      # Join the individual sequence lines into one single sequence and store in a scalar variable.
      my $sequence = join('',@seqLines); # Concatenates all elements of the @seqLines array into a single string.
      my $curr_gene = "tmp_gene";
      open TMP, ">$curr_gene" or die "Could not open $curr_gene: $!\n";
      print TMP ">$defLine\n", $sequence, "\n";
      my ($chr) = $defLine =~ m/(.*?):/;
      my $chr_file = $chr_path . $chr . ".fa";
      system("exonerate $curr_gene $chr_file -m affine:local --percent 90 --showalignment false --showvulgar false --showtargetgff true >> GENES_exonerate_output");
   }

   $/ = "\n";
   close TMP;
}


SingleFasta($genes);




# 		$cmd = "grep $curr_chr $vcf_file > tmp.vcf";
# 		system($cmd);


# exonerate gene3.5.fa CM3.6_contig30898.fa -m affine:local --showalignment false --showvulgar false --showtargetgff true

