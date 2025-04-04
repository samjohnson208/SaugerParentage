#!/usr/bin/perl
## run_first_filter_noMAF.pl by SPJ and MPR 020425
# PURPOSE: go through first filtering steps
## USAGE: perl run_first_filter_noMAF.pl [rehead.VCF]



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

my $vcf_path = "/project/ysctrout/hatchsauger/sam_sai_svit/";     ## directory containing raw vcf
my $account = "ysctrout";     ## partition
my $time = "0-00:30:00";     ## max time allowed for analysis
my $mem = "5G";             ## memory


my @misses = ('9', '8', '7', '6', '5', '4');


##########################################################
## actual work (shouldn't need to change anything below)
##########################################################

## directory for write output files
unless(-e 'first_filter_out_noMAF'){
  mkdir 'first_filter_out_noMAF', 0755 or die "Failed to make first_filter_out_noMAF directory\n";
}

my $file = shift @ARGV;


  foreach my $miss (@misses){
    my @slurmdirectives = "#!/bin/bash";
    push @slurmdirectives, "#SBATCH --job-name=filter";
    push @slurmdirectives, "#SBATCH --nodes=1";
    push @slurmdirectives, "#SBATCH --ntasks=1";
    push @slurmdirectives, "#SBATCH --cpus-per-task=1";
    push @slurmdirectives, "#SBATCH --account=$account";
    push @slurmdirectives, "#SBATCH --time=$time";
    push @slurmdirectives, "#SBATCH --mem=$mem";
    push @slurmdirectives, "#SBATCH -o "."$vcf_path"."first_filter_out_noMAF/stdout_miss"."$miss";
    push @slurmdirectives, "module load arcc/1.0 gcc/14.2.0";
    push @slurmdirectives, "cd $vcf_path";
    push @slurmdirectives, "/project/ysctrout/software/vcftools/bin/vcftools --vcf $file --out variants__miss"."$miss"." --remove-filtered-all --max-missing 0."."$miss"." --recode";

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











