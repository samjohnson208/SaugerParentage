#!/bin/bash

## slurm_cat.sh by SPJ & JPJ 110724
## PURPOSE: to cat multiple files
## USAGE: sbatch slurm_cat.sh

#SBATCH --job-name=cat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=0-4:00:00
#SBATCH --mem=1G

cd /project/ysctrout/hatchsauger/1Saug/rawreads
cat parsed_* > all_parsed.fastq

