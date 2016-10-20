#!/usr/bin/perl
#Andreu Paytuvi
#script.pl --help

if ($ARGV[0] eq "--help") {
	print "-----------------------------------------------";
	print "\nWelcome! This script will download files from database SRA and will convert from .sra format to fastq format.\n
	--output
	Specify a path where the output of the script (fastq files) will be placed. If no output is specified, files will be
	stored at the current working directory.\n\n";
	print "Call the script: './script.pl file.txt' or './script.l file.txt --output path'\n";
	print "Note1. file.txt needs to have a specific format. An example of the format:\n
	Accession	Run(s)
	ERX010637	ERR029916
	SRX118615	SRR404313	SRR404314\n\n";
	print "Note2. This script requires fastqc and fastq-dump (SRA Toolkit) in your path.\n";
	print "-----------------------------------------------\n";
	exit;
}

open FILE, $ARGV[0] or die "Could not open the file.";
$tag1 = $ARGV[1];
$outp = $ARGV[2];

if ($outp eq "") {
	$outp = "./";
}

chdir $outp;

print "Reading input\n";
sub processing_file {
    while ($lines = <FILE>) {
        next if $. < 2;
        chomp $lines;
  	@chompline = split(/\t/, $lines);
  	$accession = $chompline[0];
  	print $accession."\n";
	@runs = @chompline[1..$#chompline];
        $hash{$accession} = [ @runs ];
    }
    close FILE;
}

sub downloading_and_converting {
	print "Downloading\n";
	for $key (keys %hash) {
		mkdir $key;
		chdir $key;
		for $value ( @{$hash{$key}} ) {
			print $value."\n";
			$var1 = substr($value, 0, 3);
			$var2 = substr($value, 0, 6);
			$var3 = "wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/".$var1."/".$var2."/".$value."/".$value.".sra";
			$var4 = "fastq-dump --split-files ".$value.".sra";
			system($var3);
			print "Converting .sra file to .fastq file(s)...\n";
			system($var4);
			$var5 = "rm ".$value.".sra";
			#system($var5);
			print "Conversion performed!\n";
		}
		chdir "../";
	}
}

sub quality_analysis_fastqc {
	opendir(directory, "./");
	@dir = readdir(directory);
	closedir(directory);
	foreach $dir (@dir) {
		next if $dir =~ /\./;
		chdir $dir;
		opendir(files, "./");
		@files = readdir(files);
		closedir(files);
		mkdir fastqc_results;
		foreach $file (@files) {
			next if $file =~ /^\./;
			next if $file =~ /fastqc_results/;
			$var6 = "fastqc -o ./fastqc_results ".$file;
			system($var6);
		}
		chdir "../";
	}
}

processing_file();
downloading_and_converting();
#quality_analysis_fastqc();

print "\nDone!\n";
