#!/usr/bin/perl
## run_first_filter.pl by JPJ 20 i 22 modified by SPJ 021925
## PURPOSE: go through first filtering steps
## USAGE: perl run_first_filter_MAF1.pl [rehead.VCF]



use strict;
use warnings;
use File::Find;
use File::Path;

## usage check
if (@ARGV < 1) {
  die "Negative ghostrider, the pattern is full.\nUSAGE: perl run_first_filter.pl *bams.txt[VCF]\n";
}



##########################################################
## teton specifications to hard code
##########################################################

my $vcf_path = "/project/ysctrout/hatchsauger/sam_sai_pflav";     ## directory containing raw vcf
my $account = "ysctrout";     ## partition
my $time = "0-00:30:00";     ## max time allowed for analysis
my $mem = "5G";             ## memory

my @mafs = ('1');

##########################################################
## actual work (shouldn't need to change anything below)
##########################################################

## directory for write output files
unless(-e 'first_filter_out_MAF1'){
  mkdir 'irst_filter_out_MAF1', 0755 or die "Failed to make first_filter_out_MAF1 directory\n";
}

my $file = shift @ARGV;

foreach my $maf (@mafs){
    my @slurmdirectives = "#!/bin/bash";
    push @slurmdirectives, "#SBATCH --job-name=filter";
    push @slurmdirectives, "#SBATCH --nodes=1";
    push @slurmdirectives, "#SBATCH --ntasks=1";
    push @slurmdirectives, "#SBATCH --cpus-per-task=1";
    push @slurmdirectives, "#SBATCH --account=$account";
    push @slurmdirectives, "#SBATCH --time=$time";
    push @slurmdirectives, "#SBATCH --mem=$mem";
    push @slurmdirectives, "#SBATCH -o "."$vcf_path"."first_filter_out_MAF1/stdout_maf"."$maf";

    push @slurmdirectives, "module load arcc/1.0 gcc/14.2.0";
    push @slurmdirectives, "cd $vcf_path";
    push @slurmdirectives, "/project/ysctrout/software/vcftools/bin/vcftools --vcf $file --out variants_maf"."$maf"." --remove-filtered-all --maf 0.0"."$maf"." --recode";
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











