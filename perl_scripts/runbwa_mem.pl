#!/usr/bin/perl

# run bwa mem on fastq files passed from the command line
## USAGE: perl runbwa.pl /fullpath/to/files/*fastq

use warnings;

unless(@ARGV){
    die "runbwa.pl /fullpath/to/files/*fastq";
}

unless(-e 'sam_sai_pflav_mem'){
    mkdir 'sam_sai_pflav_mem', 0755 or die "Failed to make sam_sai_pflav directory";

foreach $fastq (@ARGV){
        $fastq =~ m/([-\w\.]+)\.fastq$/; ## fragile regexp to catch base name of individuals
        $id = $1;
        print "Mapping reads for $id\n";
        system "bwa mem -t 16 /project/ysctrout/reference_genomes/Perca_flavescens/Perca_flavescens.fasta $fastq >  aln_"."$id"."\.sam";
}


