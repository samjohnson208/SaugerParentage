#!/bin/bash

## slurm_sucker_sam2bam.sh by JPJ 7 vi 23
## PURPOSE: to use sam2bamV1.3.pl
## USAGE: sbatch slurm_sucker_sam2bam.sh

#SBATCH --job-name=sam2bam
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --account=ysctrout
#SBATCH --time=0-10:00:00
#SBATCH --mem=5G
#SBATCH -o /project/ysctrout/jjahner/suckers/slurms/std/stdout_sam2bam
#SBATCH -e /project/ysctrout/jjahner/suckers/slurms/std/stderr_sam2bam

module load arcc/1.0 gcc/12.2.0 samtools/1.16.1

cd /project/ysctrout/jjahner/suckers/sam_sai/
perl /home/jjahner/perl_scripts/sam2bamV1.3.pl *.sam
