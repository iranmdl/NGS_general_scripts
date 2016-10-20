#!/usr/bin/perl
use strict;
use warnings;
use Path::Class qw/file/;

my $file = file(shift @ARGV);

sub GenomeSize {
   my ( $genome ) = @_;
   open ( GENOME, $genome ) or die "Could not open file $genome: $!";
   # change the Input Record Separator to ">" which allows to take in a full, multiline
   # FASTA record at once.
   $/ = ">";
   my $junk = <GENOME>; # Discard the ">" at the begining of the file
   # Now read through your input file one sequence record at a time. Each input record will be a
   # multiline FASTA entry.
   my $c = 0;
   my $met = 0;
   my $stop = 0;
   my $ok = 0;
   open NOMET, ">>Melon.VALID.nomet.pep" or die "$!\n";
   open MET, ">>Melon.VALID.met.pep" or die "$!\n";
   while ( my $record = <GENOME> ) {
      chomp $record; 
      my $check;
      my ($defLine, @seqLines) = split /\n/, $record;
      # Join the individual sequence lines into one single sequence and store in a scalar variable.
      my $sequence = join('',@seqLines); # Concatenates all elements of the @seqLines array into a single string.
      my $ID = (split /\|/, $defLine)[1];
      my $first_nt = substr($sequence, 0, 1);
      my $last_nt = chop $sequence;
      $met++ if ($first_nt eq 'M');
      $stop++ if ($last_nt eq '*');
      if ( $first_nt eq 'M' && $last_nt eq '*' ) {
         $check = 'OK';
         $ok++;
         print MET ">$defLine\n$sequence$last_nt\n";
      } else {
         $check = 'ERROR';
	 if ($first_nt eq 'M' && $last_nt ne '*') {
		print NOMET ">$defLine---((VALID:met,NO_stop))\n$sequence$last_nt\n";
	 } elsif ($first_nt ne 'M' && $last_nt eq '*') {
		print NOMET ">$defLine---((VALID:NO_met,stop))\n$sequence$last_nt\n";
	 } elsif ($first_nt ne 'M' && $last_nt ne '*') {
		print NOMET ">$defLine---((VALID:NO_met,NO_stop))\n$sequence$last_nt\n";
	 }
      }
      
      open CHR_SZ, ">>numberOfPepWithM.txt" or die "Could not open numberOfPepWithM.txt: $!\n";
      print CHR_SZ "$ID\t$first_nt\t$last_nt\t$check\n";
      $c++;

   }
   print "Total number of proteins: $c\n";
   print "Total number of proteins with M: $met\n";
   print "Total number of proteins with stop: $stop\n";
   print "Total number of proteins OK: $ok\n";
   $/ = "\n";
   close (NOMET);
   close (MET);
   close ( GENOME );
}

GenomeSize($file);
