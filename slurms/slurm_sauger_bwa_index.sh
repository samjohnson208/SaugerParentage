#!/bin/bash

## slurm_sucker_bwa_index.sh by JPJ 29 v 23
## PURPOSE: bwa index
## USAGE: sbatch slurm_sucker_bwa_index.sh

#SBATCH --job-name=index
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=0-12:00:00
#SBATCH --mem=16G
#SBATCH -o /project/ysctrout/jjahner/suckers/slurms/std/stdout_bwa_index
#SBATCH -e /project/ysctrout/jjahner/suckers/slurms/std/stderr_bwa_index

module load arcc/1.0
module load gcc/12.2.0
module load bwa/0.7.17

cd /project/ysctrout/reference_genomes/c_latipinnis/
bwa index -p sucker -a bwtsw sucker_genome.fna

