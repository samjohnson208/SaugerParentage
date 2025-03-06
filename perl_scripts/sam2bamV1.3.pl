#!/usr/bin/perl
#

## Time-stamp: <Thursday, 25 February 2016, 10:55 MST -- monia>

# This script is a wrapper around samtools for batch conversion,
# sorting and indexing of sam file. This version uses fork to split
# the job over multiple processors. Adjust $ncpu for the desired
# number of cpus.

#
# Usage: perl sam2bam.pl sam_sai/*sam
#
use warnings;

if(@ARGV){
    @files = @ARGV;
}
else{
    die "Usage: sam2bam.pl sam_sai_pflav_mem/*sam";
}

my $ncpu = 4;

my $ctr = 0;
foreach (@files){
    push @{$cpu[$ctr]}, $_;
    $ctr++;
    unless($ctr % $ncpu){
	$ctr = 0;
    }
}

print "Analyzing ", scalar @files, " files\n";

foreach $i (0.. $#cpu){
    $pid = fork(); ## generate child process
    if($pid){ ## parent because fork returned non-zero
	push (@children, $pid);
    }
    elsif($pid == 0){ ## this is a child
	## walk through array of files
	foreach $j (0..$#{$cpu[$i]}){
	    $sam = $cpu[$i][$j];
	    $base = $sam;
	    $base =~ s/sam$//;
	    print "CPU ", $i+1, ", job ", $j+1, " of ";
	    print $#{$cpu[$i]}+1, ": Compressing, sorting and indexing $base\n";
	    system "samtools view -b -S -o $base"."bam $sam\n";
	    #system "samtools sort $base"."bam $base"."sorted\n"; #new in version 1.X
	    system "samtools sort $base"."bam -o $base"."sorted.bam\n";
	    system "samtools index $base"."sorted.bam\n";
	}
	exit; # end of child job
    }	
}
foreach (@children){
    waitpid($_, 0);
}
