#!/bin/bash

## slurm_sauger_split.sh by SPJ & JPJ 110724
## PURPOSE: to split fastqs by sample id
## USAGE: sbatch slurm_sauger_split.sh

#SBATCH --job-name=split
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=1-00:00:00
#SBATCH --mem=1G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_split
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_split


module load arcc/1.0 gcc/12.2.0 perl/5.34.1

cd /project/ysctrout/hatchsauger/1Saug/rawreads

perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/splitFastq_universal_regex.pl sauger_ids.txt all_parsed.fastq
