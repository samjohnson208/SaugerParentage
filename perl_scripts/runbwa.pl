#!/usr/bin/perl

## run bwa aln and samse on fastq files passed from the commandline
## USAGE: perl runbwa.pl /fullpath/to/files/*fastq

use warnings;

unless(@ARGV){
    die "runbwa.pl /fullpath/to/files/*fastq";
}

unless(-e 'sam_sai_pfluv'){
    mkdir 'sam_sai_pfluv', 0755 or die "Failed to make sam_sai_pfluv directory";
}
foreach $fastq (@ARGV){
	$fastq =~ m/([-\w\.]+)\.fastq$/; ## fragile regexp to catch base name of individuals
	$id = $1;
	print "Mapping reads for $id\n";
	system "bwa aln -n 4 -l 20 -k 2 -t 32 -R 20 -q 10 -f sam_sai_pfluv/aln_"."$id".".sai /project/ysctrout/reference_genomes/Perca_fluviatilis/europerch $fastq\n";
	## grep -c 'processor' /proc/cpuinfo  
	## tells you the number of processors on a system
	system "bwa samse -n 1 -r \'\@RG\\tID:LOC_"."$id"."\' -f sam_sai_pfluv/aln_"."$id".".sam /project/ysctrout/reference_genomes/Perca_fluviatilis/europerch sam_sai_pfluv/aln_"."$id".".sai $fastq\n";
}
