#!/usr/bin/perl
use warnings;
use strict;

my $file = shift;

sub SingleFasta {
   my ( $file ) = shift;
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
      open CHR_SZ, ">", $defLine . ".fa" or die "Could not open tmp: $!\n";
      print CHR_SZ ">$defLine\n", $sequence, "\n";
      
   }

   $/ = "\n";
   close FILE;
}


SingleFasta($file);