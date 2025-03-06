#!/bin/bash

## slurm_sauger_sam2bam.sh by JPJ and SPJ 012325
## PURPOSE: to use sam2bamV1.3.pl
## USAGE: sbatch slurm_sauger_sam2bam.sh

#SBATCH --job-name=sam2bam
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --account=ysctrout
#SBATCH --time=0-10:00:00
#SBATCH --mem=5G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_sam2bam
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_sam2bam

module load arcc/1.0 gcc/14.2.0 samtools/1.20

cd /project/ysctrout/hatchsauger/sam_sai_pflav_mem/
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/sam2bamV1.3.pl *.sam

