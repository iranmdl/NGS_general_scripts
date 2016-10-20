#!/usr/bin/perl -w

# Irantzu Anzar
# 20/03/2015

use 5.010;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use POSIX qw(strftime);
use List::Util 'first';


# Key as the normal (wt) genotype and values contain the mutant genotypes to be discarded
# for each wt genotype
my %Genotypes = (
    "0/0"  => [ "./.", "0/0" ],
    "0/1"  => [ "./.", "0/0", "0/1" ],
    "1/1"  => [ "./.", "0/0", "0/1", "1/1" ],
    "0/2"  => [ "./.", "0/0", "0/2" ],
    "1/2"  => [ "./.", "0/0", "0/1", "0/2", "1/2" ],
    "2/2"  => [ "./.", "0/0", "0/2" ],
    "0/3"  => [ "./.", "0/0", "0/3" ],
    "1/3"  => [ "./.", "0/0", "0/3", "0/1", "1/3" ],
    "2/3"  => [ "./.", "0/0", "0/3", "0/2", "2/3" ],
    "3/3"  => [ "./.", "0/0", "0/3" ],
);


my ($vcf, $normal, $mut, $help, $man);

GetOptions (
         'vcf|v=s' => \$vcf,
         'normal|n=s' => \$normal,
         'mutant|m=s' => \$mut,
         'help|h' => \$help,
         'man' => \$man,
) or pod2usage(-verbose => 1) && exit;


pod2usage(-verbose => 1) && exit if defined $help;
pod2usage(-verbose => 2) if defined $man;
pod2usage("$0: No files given.") if ((!defined $vcf) && (-t STDIN));

#   Check for required variables.
if ( !($normal) or !($mut) ) { 
    pod2usage({EXITSTATUS => 2, VERBOSE => 1, MESSAGE => "\n$0: Missing arguments.\n",});
}

# Output fixed to vcf folder
my $outdir = dirname($vcf);
# Output file name and path
my $bs = basename($vcf);
$bs =~ s{\.[^.]+$}{};   # removes extension
my $out = "$outdir/$bs.$mut.VS.$normal.vcf";
# Remove files if exist
unlink ($out) if (-e $out);



#----------------------------------------------------------------------
#---------------------  Start variant annotation ----------------------

print "\n========= Variant Annotation Script =========\n";
printf '=' x45;

print
	"\n>\tInput file: $vcf\n",
	">\tOutput file: $out\n",
	">\tComparison: $mut VS $normal\n",
	">\tAnalysis started - ", strftime("[%H:%M:%S]", localtime),
    "\n";




my ( @record, $index_n, $index_m, $gt_n, $gt_m );
my ( $c1, $c2 )  = ( 0, 0 );

open ( my $vcf_fh,  "<$vcf" ) or die "Unable to open file $vcf: $!";
open ( my $out_fh, ">>$out" ) or die "Unable to open $out: $!\n";

while ( <$vcf_fh> ) {
	@record = split;
	chomp(@record);
	if ( $_ =~ /^\##/ ) {
		print $out_fh join("\t", @record), "\n";
	}
	elsif ( $_ =~ /^\#CHROM/ ) {
		( $index_n ) = first { $record[$_] eq $normal } 0 .. $#record;
		( $index_m ) =  first { $record[$_] eq $mut } 0 .. $#record;
		print $out_fh join("\t", @record), "\n";
		# print $out_fh join("\t", @record[0..8]), "\t$record[$index_n]\t$record[$index_m]\n";
	}
	elsif ( $_ ne /^\#/ || $_ ne /^\s*$/ ) {
		if ( ! defined $index_n ) {
			print "Error. Input vcf file without header or the provided sample names are not equal to 
					the names present in vcf header\n";
			exit;
		}
		print "\r\t   ... Processed $c1 entries ...";
		# my @fields = @record[0..8];
		$gt_n = (split /:/, $record[$index_n])[0];
		$gt_m = (split /:/, $record[$index_m])[0];
		if ( exists $Genotypes{$gt_n} ) {
			if ( !grep( /^$gt_m$/, @{$Genotypes{$gt_n}} ) ) {
				# print $out_fh join("\t", @fields),"\t$record[$index_n]\t$record[$index_m]\n";
				print $out_fh join("\t", @record),"\n";
				$c2++;
			}
		}
		$c1++;
	}
}

close($vcf_fh);
close($out_fh);


print
	"\n>\tAnalysis finished - ", strftime("[%H:%M:%S]", localtime),
	"\n>\tTotal records processed: $c1\n",
	">\tTotal records kept: $c2\n",
    "\n";