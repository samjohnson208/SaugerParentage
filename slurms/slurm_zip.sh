#!/bin/bash

## slurm_zip.sh by JPJ 27 i 23
## PURPOSE: to zip or unzip big files
## USAGE: sbatch slurm_zip.sh

#SBATCH --job-name=zip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

gzip /project/ysctrout/hatchsauger/1Saug/rawreads/1SaugEvens.fastq
gzip /project/ysctrout/hatchsauger/1Saug/rawreads/1SaugOdds.fastq


