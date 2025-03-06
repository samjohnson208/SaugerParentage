#!/usr/bin/perl

# run bwa mem on fastq files passed from the command line
## USAGE: perl runbwa.pl /fullpath/to/files/*fastq

use warnings;

unless(@ARGV){
    die "runbwa.pl /fullpath/to/files/*fastq";
}

unless(-e 'sam_sai_pflav_mem'){
    mkdir 'sam_sai_pflav_mem', 0755 or die "Failed to make sam_sai_pflav_mem directory";
}

my $fastq;
foreach $fastq (@ARGV){

    $fastq =~ m/([A-Z]+[0-9]+_[0-9]+)\.fastq/;
    #print $fastq;
    my $id = $1;
    my $sam = "aln_"."$id"."\.sam";
    my $bam = "aln_"."$id"."\.bam";
    my $sorted = "aln_"."$id"."\.sorted\.bam";

    if(-e "/project/ysctrout/hatchsauger/sam_sai_pflav_mem/$sorted"){
	print "$sorted already exists\n";
    }else{
    print "Mapping reads for $id\n";

push @jobarray, "bwa mem -t 16 /project/ysctrout/reference_genomes/Perca_flavescens/Perca_flavescens.fasta $fastq >  aln_$id.sam \n 

	echo \"Converting sam to bam for $id\n\" \n

	samtools view -b -S -o $bam $sam \n

	echo \"Sorting and indexing bam files for $id\n\" \n
	samtools sort $bam -o $sorted \n
	samtools index $sorted
 \n";
    }
}


