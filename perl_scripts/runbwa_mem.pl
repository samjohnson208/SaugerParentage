#!/usr/bin/perl

# run bwa mem on fastq files passed from the command line. run from an empty sam_sai_species_mem directory
## USAGE: perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/runbwa_mem.pl /project/ysctrout/hatchsauger/1Saug/rawfastqs/*.fastq
## USAGE: perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/runbwa_mem.pl /project/ysctrout/hatchsauger/1Saug/rawfastqs/SAR_16_6550.fastq
## USAGE: perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/runbwa_mem.pl /project/ysctrout/hatchsauger/1Saug/contam_cleaned/ind_files/SAR_21_6100.fastq


use warnings;

my @jobarray = "#!/bin/bash";

unless(@ARGV){
    die "runbwa.pl /fullpath/to/files/*fastq";
}

my $fastq;
foreach $fastq (@ARGV){

    # $fastq =~ m/([A-Z]+_[0-9]+_[0-9]+)\.fastq/;
    # my $id = $1;
    
    # THIS CHUNK ABOVE THAT I USUALLY USE DID NOT WORK ON CONTAM TO FASTP FASTQS
    # because instead of just being SAR_YY_XXXX.fastq, they're trim_read_SAR_YY_XXXX.fastq
    # after being output from fastp. This new chunk (27 - 35) should fix this.
    
    # Remove path and extension
	(my $file = $fastq) =~ s!.*/!!;
	$file =~ s/\.fastq$//;

	# Extract SAR_Year_SampleNumber
	my ($id) = $file =~ /(SAR_\d+_\d+)/;
    my $sam = "mem_svit_contam_fastp_"."$id"."\.sam";
    my $bam = "mem_svit_contam_fastp_"."$id"."\.bam";
    my $sorted = "mem_svit_contam_fastp_"."$id"."\.sorted\.bam";


push @jobarray, "#SBATCH --account=ysctrout";
push @jobarray, "#SBATCH --job-name=bwa_mem";
push @jobarray, "#SBATCH --time=7-00:00:00"; 
push @jobarray, "#SBATCH --nodes=1";
push @jobarray, "#SBATCH --ntasks-per-node=32"; # one core per node
push @jobarray, "#SBATCH --mem=64000";
push @jobarray, "#SBATCH --mail-type=END";
push @jobarray, 'module load arcc/1.0 gcc/14.2.0 bwa/0.7.17 samtools/1.20'; 

push @jobarray, "bwa mem -t 32 /project/ysctrout/reference_genomes/Sander_vitreus/walleye $fastq >  mem_svit_contam_fastp_"."$id".".sam"; 

push @jobarray, "echo \"Converting sam to bam for "."$id"."\n\"";

push @jobarray, "samtools view -b -S -o $bam $sam";

push @jobarray, "echo \"Sorting and indexing bam files for "."$id"."\n\"";
	
push @jobarray, "samtools sort $bam -o $sorted";
	
push @jobarray, "samtools index $sorted";

    }


my $slurm = join "\n", @jobarray;
# print $slurm\n;

## final job
runserialjob($slurm);

#### -------------------------------------------------------------------
sub runserialjob{
    my $slurmjob = $_[0];
    $slurmjob .= "\nexit\n";
    open SBATCH, "| sbatch 1>/dev/null" or die "Failed to fork for sbatch; $!\n";
    print SBATCH "$slurmjob";
    close(SBATCH) or die "Couldn't close SBATCH\n";
}




