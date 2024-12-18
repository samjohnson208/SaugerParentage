#!/usr/bin/perl
## run_parse_evens.pl by SPJ & JPJ 110424
## PURPOSE: parse multiple fastq files
## USAGE: perl run_parse_evens.pl even_fq_*


use strict;
use warnings;
use File::Find;
use File::Path;

## usage check
if (@ARGV < 1) {
  die "Negative ghostrider, the pattern is full.\nUSAGE: perl run_parse.pl split_fq_*\n";
}



##########################################################
## teton specifications to hard code
##########################################################

my $account = "ysctrout";     ## partition
my $time = "5-00:00:00";     ## max time allowed for analysis
my $mem = "1G";             ## memory


##########################################################
## actual work (shouldn't need to change anything below)
##########################################################



## launch new job for each input file

foreach my $file (@ARGV) {
  ## standard sbatch
  
  my @slurmdirectives = "#!/bin/bash";
  push @slurmdirectives, "#SBATCH --job-name=parse";
  push @slurmdirectives, "#SBATCH --nodes=1";
  push @slurmdirectives, "#SBATCH --ntasks=1";
  push @slurmdirectives, "#SBATCH --cpus-per-task=1";
  push @slurmdirectives, "#SBATCH --account=$account";
  push @slurmdirectives, "#SBATCH --time=$time";
  push @slurmdirectives, "#SBATCH --mem=$mem";

  
  ## actual work
  
  push @slurmdirectives, "module load arcc/1.0 gcc/14.2.0";
  push @slurmdirectives, "cd /project/ysctrout/hatchsauger/1Saug/rawreads";
  push @slurmdirectives, "~/miniconda3/bin/perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/parse_barcodes768.pl 1SaugEvens_Demux.csv $file LH0";

  
  ## join slurmdirectives and print
  
  my $slurm = join "\n", @slurmdirectives;
  runserialjob($slurm);
  #print "$slurm\n";
}

print "Yippee-ki-yay, motherfucker!\n";


## subroutine that initializes the slurm job

sub runserialjob{
  my $slurmjob = $_[0];
  $slurmjob .= "\nexit\n";
  open SBATCH, "| sbatch 1>/dev/null" or die "Failed to fork for sbatch; $!\n";
  print SBATCH "$slurmjob";
  close(SBATCH) or die "Couldn't close SBATCH\n";
}


