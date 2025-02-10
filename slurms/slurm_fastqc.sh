#!/bin/sh

## slurm_fastqc.sh by SPJ on 021024, modified by SPJ 121924
## PURPOSE: to check the quality of raw reads for each library (large, unzipped .fastq files)
## USAGE: sbatch slurm_fastqc.sh

#SBATCH --job-name=fastqc
#SBATCH --account=ysctrout
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15G
#SBATCH --mail-type=END


module load arcc/1.0
module load fastqc/0.12.1

cd /project/ysctrout/hatchsauger/1Saug/rawreads

mkdir -p "out_fastqc"

fastqc *.fastq -o "$out_fastqc"

