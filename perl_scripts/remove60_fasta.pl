#!/usr/bin/perl

##This script removes line endings after 60 characters in fasta file

$in = shift(@ARGV);

open (IN, "< $in") or die "couldnt open IN\n";
open (OUT, "> no60_$in") or die "couldnt open OUT\n";

while(<IN>){
    if ($_ =~ m/>/){
	print OUT "\n"."$_";
    }
    elsif ($_ =~ m/[ATCG]/){
	$_ =~ s/\n//;
	print OUT "$_";
    }

}

close(IN);
close(OUT);
