#!/bin/bash

## slurm_cat.sh by JPJ 28 i 23
## PURPOSE: to cat multiple files
## USAGE: sbatch slurm_cat.sh

#SBATCH --job-name=cat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=evolgen
#SBATCH --time=0-1:00:00
#SBATCH --mem=1G

cd /project/evolgen/jjahner/sab/rawreads/
cat parsed_split_fq_* > all_parsed.fastq
