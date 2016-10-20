#!/usr/bin/perl
use strict;
use warnings;
use LWP::Simple;

my $file = "new_IDs_tab.txt";
#my $outfile = "KEGG_parsed_out.txt"


open (FH, $file) or die "Could not open file $file: $!";
open (OUTPUT, ">>KEGG_parsed_out_newIDs.txt" ) or die "Could not open: $!\n";

while (<FH>) {
	chomp;
    my @lines = split /\t/;
    my $ID = $lines[0];
    print "Processing $ID ...\n";

	# Web pages
	my $firsturl = "http://www.genome.jp/dbget-bin/get_linkdb?-t+uniprot+";
	my $secondurl = "http://www.genome.jp/dbget-bin/www_bget?uniprot:";
	
	# Web page with current ID
	my $first_url_ID = $firsturl . $ID;
	  
	my $first_WEB = get $first_url_ID;
		die "Couldn't get $first_url_ID" unless defined $first_WEB;

	#print "$content";
	my $firstID = $1 if ($first_WEB =~ m/\?uniprot\:(.*?)">/);
	
	if ( ! defined $firstID ) {
    	print OUTPUT join("\t", @lines),"\tNone\n";
    	next;
    }
	
	# print "$firstID\n";
	my $sec_url_ID = $secondurl . $firstID;

	my $second_WEB = get $sec_url_ID;
		die "Couldn't get $sec_url_ID" unless defined $second_WEB;

	my $secID = $1 if ($second_WEB =~ m/EnsemblPlants; (.*?);/);
	if ( ! defined $secID ) {
    	print OUTPUT join("\t", @lines),"\tNone\n";
    	next;
    }

	# print "$secID\n";
	print OUTPUT join("\t", @lines),"\t$secID\n";

}

close (FH);
close (OUTPUT);








# http://www.genome.jp/dbget-bin/get_linkdb?-t+uniprot+sly:100037510

# http://www.genome.jp/dbget-bin/www_bget?uniprot:A1IKU3