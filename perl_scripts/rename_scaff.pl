#!/usr/bin/perl

##script to rename identifier lines in genome assembly .fasta file so that they all have consistent names. Assumes no line endings within sequences (no60).

#usage: perl rename_scaff.pl assembly.fasta



my $in = $ARGV[0];

open(IN, "< $in") or die "couldnt open IN\n";
open(OUT, "> renamed_".$in.".txt") or die "couldnt open OUT\n";
open(OUT2, "> rename_scaff_key.txt") or die "couldnt open OUT2\n";

my $ctr=0;


while(<IN>){
    chomp;
    if(/^>/){
          $ctr++;
	  print OUT ">scaffold_".$ctr."\n";
          print OUT2 "$_"." scaffold_".$ctr."\n";
    }
    else{
      print OUT "$_\n";
    }
}

close(IN);
close(OUT);
close(OUT2);

