#!/bin/bash

## slurm_zip.sh by JPJ 27 i 23
## PURPOSE: to zip big files
## USAGE: sbatch slurm_zip.sh

#SBATCH --job-name=zip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=1-00:00:00
#SBATCH --mem=1G


gunzip /project/ysctrout/hatchsauger/1Saug/rawreads/*.fastq.gz

