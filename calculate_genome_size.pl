#!/usr/bin/perl

if( $ARGV[0] eq '-h' || $ARGV[0] eq '-help') {
	help();
	exit;
}


sub help { print "USAGE:\n\tGenomeSize [ FILE.fasta] [-help]\n\n";}

$/ = ">";
my $junk = <>;
my $g_size = 0;   
while ( my $record = <> ) {
   chomp $record; 
   my ($defLine, @seqLines) = split /\n/, $record; 
   my $sequence = join('',@seqLines);
   $g_size+=length($sequence);
   print "$defLine\t", length($sequence), "\n";  
}
$/ = "\n";
