#!/bin/bash

## slurm_sucker_split.sh by JPJ 31 v 23
## PURPOSE: to split fastqs by ind
## USAGE: sbatch slurm_sucker_split.sh

#SBATCH --job-name=split
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=1-00:00:00
#SBATCH --mem=1G
#SBATCH -o /project/ysctrout/jjahner/suckers/slurms/std/stdout_split
#SBATCH -e /project/ysctrout/jjahner/suckers/slurms/std/stderr_split


module load arcc/1.0 gcc/12.2.0 perl/5.34.1

cd /project/ysctrout/jjahner/suckers/rawreads/

perl /home/jjahner/perl_scripts/splitFastq_universal_regex.pl sucker_ids.txt all_parsed.fastq
