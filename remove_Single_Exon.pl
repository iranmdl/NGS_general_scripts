#!/usr/bin/perl -w

# Irantzu Anzar
# 04/08/2015


use 5.010;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Data::Dumper;
use List::Util qw(first);


my ($gtf, $help, $man);

GetOptions (
    	'gtf|g=s' => \$gtf,
    	'help|h' => \$help,
    	'man' => \$man,
) or pod2usage({VERBOSE => 0, MESSAGE => "Invalid arguments.\n"});
   
pod2usage(VERBOSE => 1) && exit if defined $help;
pod2usage(VERBOSE => 2) if defined $man;
#   Check for required variables.
pod2usage({EXITSTATUS => 2, VERBOSE => 1, MESSAGE => "\n$0: Missing arguments.\n",}) if !($gtf);
pod2usage({EXITSTATUS => 2, VERBOSE => 1, MESSAGE => "\n$0: $gtf does not exist.\n",}) unless ( -e $gtf );


######### Output ###########
# Output fixed to annotation folder
my $outdir = dirname($gtf);
# Output file name and path
my $bs = basename($gtf);
$bs =~ s{\.[^.]+$}{}; 			# removes extension
my $out_gtf = "$outdir/$bs.rmSE.gtf";
my $log_f = "$outdir/Remove_SE_script.log";
my $transcript_list = "$outdir/Kept_transcript_IDs.log";
# Remove files if exist
unlink ($out_gtf) if (-e $out_gtf);
unlink ($log_f) if (-e $log_f);
unlink ($transcript_list) if (-e $transcript_list);

#----------------------------------------------------------------------
#----------------------------  Subroutines ----------------------------
sub gtf_extraction {
	my ( $gtf, $log_f ) = @_;
	print "1. Reading gtf file ...\n";
	open my $gtf_f, '<', $gtf or die "Can't read from $gtf: $!";
	open my $log, ">", "$log_f" or die "$!\n";
	my %gtf_info;
	while (<$gtf_f>){
		chomp;
    	next if ($_ =~ /^\#/ || $_ =~ /^\s*$/); 
		my @field = split /\t/;
		my ( $gene, $transcript ) = (split /;/, $field[8])[0,1];
    	$gene =~ s/gene_id "//g;
    	$transcript =~ s/ transcript_id "//g;
    	foreach ( $gene, $transcript ) { $_ =~ s/"//g; };
    	if ( $field[6] eq '+' ) {
    		$transcript .= '_PLUS';
    	}
    	elsif ( $field[6] eq '-' ) {
    		$transcript .= '_NEG';	
    	}
		if ( $field[2] eq 'exon' ) {
			push @{$gtf_info{$gene}{$transcript}}, 1;
		}
    }
    print $log Dumper \%gtf_info;
    close($log);
    close($gtf_f);
    return(\%gtf_info);
}

sub gtf_parsing {
	my ( $gtf_info ) = @_;
	print "2. Processing gtf file ...\n";
	my %gtf_info = %$gtf_info;
	my @kept_transcripts;

	foreach my $gene ( sort keys %gtf_info ) {
		my $transcript_id;
		my $tr_strand;
		my %to_evaluate;
		my @strands;
		foreach my $transc ( sort keys %{$gtf_info{$gene}} ) {
			if ( index($transc, '_PLUS') != -1 ) {
	    		( $transcript_id = $transc ) =~ s/_PLUS//g;
	    		$tr_strand = 1;
			}
			else {
				( $transcript_id = $transc ) =~ s/_NEG//g;
				$tr_strand = 2;
			}

			if ( scalar(keys %{$gtf_info{$gene}}) == 1 ) {      # Only one isoform: Keep transcript
				push @kept_transcripts, $transcript_id;
			}
			else {
				my @exons = @{$gtf_info{$gene}{$transc}};
				if ( scalar@exons > 1 ) {
					push @kept_transcripts, $transcript_id;
					push @strands, $tr_strand;
				}
				else {
					$to_evaluate{$transcript_id} = $tr_strand;
				}
				foreach my $evaluated_transcript ( sort keys %to_evaluate ) {
					if ( first { $_ eq $to_evaluate{$evaluated_transcript} } @strands ) {
						next;
					}
					else {
						push @kept_transcripts, $transcript_id;
					}
				}
			}
		}
	}
	return( @kept_transcripts );
}

sub print_output {
	my ( $out_gtf, $transcript_list, $gtf, @kept_transcripts ) = @_;
	print "3. Printing output file ...\n";
	open LIST, ">", "$transcript_list" or die "$!\n";
	print LIST join( "\n", @kept_transcripts ), "\n";
	system ("grep -w -F -f $transcript_list $gtf > $out_gtf");
	close(LIST);
}  


my ( $gtf_info ) = gtf_extraction( $gtf, $log_f );
my ( @kept_transcripts ) = gtf_parsing($gtf_info);
print_output( $out_gtf, $transcript_list, $gtf, @kept_transcripts );

# print join( "\n", @kept_transcripts ), "\n";

