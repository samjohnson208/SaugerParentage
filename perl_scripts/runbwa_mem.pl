#!/usr/bin/perl

# run bwa mem on fastq files passed from the command line
## USAGE: perl runbwa.pl /fullpath/to/files/*fastq

use warnings;

my @jobarray = ();

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

push @jobarray, "#SBATCH --account=ysctrout";
push @jobarray, "#SBATCH --job-name=bwa_mem";
push @jobarray, "#SBATCH --time=2:00:00"; 
push @jobarray, "#SBATCH --nodes=1";
push @jobarray, "#SBATCH --ntasks-per-node=16"; # one core per node
push @jobarray, "#SBATCH --mem=64000"; 
push @jobarray, "module load arcc/1.0 gcc/14.2.0 bwa/0.7.17"; 

push @jobarray, "bwa mem -t 16 /project/ysctrout/reference_genomes/Perca_flavescens/Perca_flavescens.fasta $fastq >  aln_"."$id".".sam \n 

	echo \"Converting sam to bam for "."$id"."\n\" \n

	samtools view -b -S -o $bam $sam \n

	echo \"Sorting and indexing bam files for "."$id"."\n\" \n
	samtools sort $bam -o $sorted \n
	samtools index $sorted
 \n";
    }
}

my $slurm = join "\n", @jobarray;
print $slurm\n;

## final job
# runserialjob($slurm);

#### -------------------------------------------------------------------
sub runserialjob{
    my $j = $_[0];
    my $slurmjob = '';
    $slurmjob .= $pbsconf;
    $slurmjob .= $prolog;
    $slurmjob .= $jobarray[$j];
    $slurmjob .= $epilog;
    $slurmjob .= "exit\n";
    open SBATCH, "| sbatch 1>/dev/null" or die "Failed to fork for sbatch; $!";
    print SBATCH "$slurmjob";
    close(SBATCH) or die "Couldn't close SBATCH";
}






