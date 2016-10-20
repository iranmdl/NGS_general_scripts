#!/usr/bin/perl -w

# Irantzu Anzar
# 02/07/2015


use 5.010;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Data::Dumper;
use List::MoreUtils qw(any);


my ($des, $gtf, $of, $gene_num, $help, $man);

GetOptions (
    	'described|d=s' => \$des,
    	'gtf|g=s' => \$gtf,
    	'official_name_prefix|of=s' => \$of,
    	'last_gene_number|n=i' => \$gene_num,
    	'help|h' => \$help,
    	'man' => \$man,
) or pod2usage({VERBOSE => 0, MESSAGE => "Invalid arguments.\n"});
   
pod2usage(VERBOSE => 1) && exit if defined $help;
pod2usage(VERBOSE => 2) if defined $man;
#   Check for required variables.
pod2usage({EXITSTATUS => 2, VERBOSE => 1, MESSAGE => "\n$0: Missing arguments.\n",}) if ( !($des) or !($gtf) or !($of) or !($gene_num) );
pod2usage({EXITSTATUS => 2, VERBOSE => 1, MESSAGE => "\n$0: $des or $gtf does not exist.\n",}) unless ( -e $des || -e $gtf );

######### Output ###########
# Output fixed to annotation folder
my $outdir = dirname($gtf);
# Output file name and path
my $bs1 = basename($gtf);
my $bs2 = basename($des);
foreach ( $bs1, $bs2 ) { $_ =~ s{\.[^.]+$}{}; }; 			# removes extension
my $out_gtf = "$outdir/$bs1.reformat.gtf";
my $out_des = "$outdir/$bs2.reformat.txt";
my $log_f = "$outdir/PASA_reformat.log";
# Remove files if exist
unlink ($out_gtf) if (-e $out_gtf);
unlink ($out_des) if (-e $out_des);
unlink ($log_f) if (-e $log_f);

#----------------------------------------------------------------------
#----------------------------  Subroutines ----------------------------

# sub rename_duplicates {
# 	my $hash = shift;
# 	my %reverse;
# 	open my $fh, '>>', $out_dup_log
# 		or die "$0: cannot open $out_dup_log: $!\n";
# 	print $fh "# Duplicate transcript names in input file.\nThey will be tagged with 'X1', 'X2' ...\n";
# 	while (my ($key, $value) = each %{$hash}) {
# 		push @{$reverse{$value}}, $key;
# 	}
# 	while (my ($key, $value) = each %reverse) {
# 		next unless @$value > 1;
# 		print $fh scalar(@$value), " duplicates found: @$value have the same value $key\n";
# 		my $c = 1;
# 		for my $rep ( 0 .. $#$value ) {
# 			next if $rep == 0;
# 			$hash->{ @$value[$rep] } = $key . 'X' . $c;
# 			$c++;
# 		}
# 	}
# }


sub Extract_tag {
	my ( @all_n ) = @_;

	my ( $of_c, $cuff_c, $est_c ) = ( 0, 0, 0 );
	# $of_c = 1 if ( grep /^$of/, @all_n );
	# $cuff_c = 1 if ( grep /^TCONS/, @all_n );
	# $est_c = 1 if ( ( grep !/^$of/, @all_n ) && ( grep !/^TCONS/, @all_n ) )
	
	if ( grep /^$of/, @all_n ) {
		@all_n = grep {!/^$of/}@all_n;
		$of_c = 1;
	}
	if ( grep /^TCONS/, @all_n ) {
		@all_n = grep {!/^TCONS/}@all_n;
		$cuff_c = 1;
	}
	if ( scalar@all_n > 0 ) {
		$est_c = 1;
	}

	# Check all options
	my $tag;
	if ( ($of_c == $cuff_c) && ($of_c == 1) && ($est_c == 1) ) {
		$tag = " official \"yes\"; RNASeq \"yes\"; EST \"yes\"";
	}
	elsif ( ($of_c == $cuff_c) && ($of_c == 1) && ($est_c == 0) ) {
		$tag = " official \"yes\"; RNASeq \"yes\"; EST \"no\"";
	}		
	elsif ( ($of_c == $cuff_c) && ($of_c == 0) && ($est_c == 1) ) {
		$tag = " official \"no\"; RNASeq \"no\"; EST \"yes\"";
	}
	elsif ( ($of_c == $cuff_c) && ($of_c == 0) && ($est_c == 0) ) {
		$tag = " official \"no\"; RNASeq \"no\"; EST \"no\"";# Warn! not possible.
		print "ERROR.\nTranscript without any support found. Exiting...\n" && exit;
	}
	elsif ( ($of_c != $cuff_c) && ($of_c == 1) && ($est_c == 0) ) {
		$tag = " official \"yes\"; RNASeq \"no\"; EST \"no\"";
	}
	elsif ( ($of_c != $cuff_c) && ($of_c == 0) && ($est_c == 0) ) {
		$tag = " official \"no\"; RNASeq \"yes\"; EST \"no\"";
	}
	elsif ( ($of_c != $cuff_c) && ($of_c == 1) && ($est_c == 1) ) {
		$tag = " official \"yes\"; RNASeq \"no\"; EST \"yes\"";
	}
	elsif ( ($of_c != $cuff_c) && ($of_c == 0) && ($est_c == 1) ) {
		$tag = " official \"no\"; RNASeq \"yes\"; EST \"yes\"";
	}
	return $tag;
}




#----------------------------------------------------------------------
#---------------------  Parse described text file ---------------------

# print "\n============ PASA_reformat ============\n";
# printf '=' x40;
# print "\n   > Parsing peptide fasta file:\n";

# Open described file and save into array
open my $des_fh, '<', $des
	or die "$0: cannot access $des: $!";
chomp( my @lines_des = <$des_fh> );
close $des_fh;

# Open output files
open my $out_des_fh, '>>', $out_des
	or die "$0: cannot open $out_des: $!\n";
# Log file
open my $log_fh, '>>', $log_f
	or die "$0: cannot open $log_f: $!\n";

my ( %transcript_H, %locus_H, %joined_H, @transcript_lines, @all_ids, @all_gene_ids );

for my $nr ( 0 .. $#lines_des ) {  
	next if ( $lines_des[$nr] =~ /^\s+$/ ); 									# Skip blank lines
    my @nr_arr = split /\t/, $lines_des[$nr];
    my $last = pop @nr_arr;
    my $line = join "\t", @nr_arr; 		# All line except alignment description field
    
    # Header
    if ( $lines_des[$nr] =~ /^#/ ) {
		print $out_des_fh join("\t", @nr_arr),"\tTranscript_official \tLocus_official\n";
		next;
    }

    # Capture interesting columns
    my ( $locus, $names ) = (@nr_arr)[1,3];
    # Evaluate next line
    my $next_line_locus;
    if ( $nr != $#lines_des ) {
		$next_line_locus = (split /\t/, $lines_des[$nr + 1])[1];
    }
    else {
    	$next_line_locus = '10000000000000000';
    }
	
	my @curr_ids = split( /,/, $names );
	push @transcript_lines, $line;
	push @all_ids, @curr_ids;

	if ( $next_line_locus != $locus ) {	    						# Next line corresponds to another locus, so we need to print current locus results
		my @official_names = grep { /^$of/ } @all_ids;
		my $gene_id;
		if ( scalar@official_names >= 1 ) {							# Extract official names
			$gene_id = (sort @official_names)[0];
			$gene_id =~ s{\.[^.]+$}{};								# Remove all after dot if a dot is present
			if ( grep( !/^$of/, @all_ids ) ) { 							# The locus is also supported by an EST or RNASeq ---> Add 'C' to name (corrected)
				$gene_id .= 'C';
			}
		}
		else {														# If there is not official  names, look for cufflinks names
			my @cuff_names = grep { /^TCONS/ } @all_ids;
			if ( scalar@cuff_names >= 1 ) {
				$gene_num++;
				if ( $gene_num < 100000 ) {    # Porque así si el numero más grande es 027427, nosotros al programa le damos 27427, y el imprimirá 027428 para el siguiente, y no 27428
					$gene_id = $of . '0' . $gene_num;
				}
				else {
					$gene_id = $of . $gene_num;
				}
			}
			else {
				print $log_fh "Warn! $locus is a EST prediction. The gene should be supported by cufflinks or official transcript to be processed...\n";			# This transcript is only supported by EST ---> We don't trust it
				splice( @all_ids );
				splice( @transcript_lines );
				next;
			}
    	}
    	
    	# Locus hash
    	my $locus_name = "PASA_cluster_" . $locus;
    	if ( $gene_id ~~ @all_gene_ids ) {
			print $log_fh "Houston we have a problem. Repeated gene_id $gene_id in \n";
			$gene_id .= '_bis';
		}
    	push @all_gene_ids, $gene_id;
    	$locus_H{$locus_name} = $gene_id;	
    	
    	my $transcript_num = 1;
    	my $transcript_name;
    	my $joined;
    	for my $transcript_line ( @transcript_lines ) {
			$transcript_name = $gene_id . '.' . $transcript_num;
			print $out_des_fh $transcript_line,"\t$transcript_name\t$gene_id\n";	
			my $asmbl_id = (split /\t/, $transcript_line)[2];
			$transcript_H{$asmbl_id} = $transcript_name;
			my $joined = (split /\t/, $transcript_line)[3];
			my @joined = (split /,/, $joined);
			$joined_H{$asmbl_id} = $joined;
			$transcript_num++;
		}
		splice( @all_ids );
		splice( @transcript_lines );
	}
}

close $out_des_fh;
close $log_f;
# # print Dumper( \%locus_H );
# # print Dumper( \%transcript_H );
# print Dumper( \%joined_H );

# # Check if there are duplicates in transcript names
# # If so, rename them with 'X1', 'X2' ...
# rename_duplicates(\%locus_H);
# rename_duplicates(\%transcript_H);

# Merge both hashes for easier substitution in gtf file
my %merged_H = ( %locus_H, %transcript_H );

# open my $fh, '>>', $out_dup_log
#                 or die "$0: cannot open $out_dup_log: $!\n";
# print $fh Dumper( \%merged_H );
# print $fh "\n##########\n";
# print $fh Dumper( \%joined_H );
# close($fh);


# Open gtf file
open my $gtf_fh, '<', $gtf
    or die "$0: cannot access $gtf: $!\n";
open my $out_gtf_fh, '>>', $out_gtf
	or die "$0: cannot open $out_gtf: $!\n";

# Replacement of PASA names with official names
my $re = join '|', map quotemeta, keys %merged_H;
my $join = join '|', map quotemeta, keys %joined_H;

while ( <$gtf_fh> ) {
	chomp;
	my $tag;
	if ( /\btranscript\b/ && /\b($join)\b/ ) {
		my @all_n = split /,/, $joined_H{$1};
		my $tag = Extract_tag(@all_n);
		$_ = $_ . $tag;
		$_ = $_ . "; Coming_from \"$joined_H{$1}\"";
		# $_ = $_ . "Joined:$joined_H{$1}";
	}
	s/\b($re)\b/$merged_H{$1}/g;
	print $out_gtf_fh $_, "\n";
}

close $gtf_fh;
close $out_gtf_fh;



__END__
 
=head1 NAME
 
PASA_reformat.pl - PASA gtf reformat script
 
=head1 SYNOPSIS

 perl PASA_reformat.pl [--help|-h] [--man]
                  --described|-d XXX.pasa_assemblies_described.txt FILE
                  --gtf|-g XXX.pasa_assemblies.gtf FILE
                  --official _name|-of NAME

=head1 ARGUMENTS

 --described|-d    File with assembly (clusters) information ( XXX.pasa_assemblies_described.txt) filename and path
 --gtf|-g    Annotation (XXX.pasa_assemblies.gtf) filename and path
 --official_name_prefix|-of    Official name gene prefix 
 --last_gene_number|-n    Biggest gene number in official annotation

=head1 OPTIONS

 --help        Print a brief help message and exits.
 --man         Prints the manual page and exits.

=head1 DESCRIPTION
 
 This program will read the given assembly description file from PASA
 output and will replace PASA cluster names with official  transcript names or
 cufflinks names for de-novo predicted transcripts. The output is saved
 in the same directory as the described txt file with input file name plus
 'reformat' tag.
 ATTENTION!
 ===> Sort described file by locus id before running PASA_reformat.pl <===
 awk 'NR == 1; NR > 1 {print $0 | "sort -k2,2n"}' XXX.described.txt > XXX.described_sort.txt
 
=cut

