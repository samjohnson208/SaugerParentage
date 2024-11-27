#!/bin/bash

## slurm_woodcock_runbwa.sh by JPJ 5 vi 23
## PURPOSE: to use runbwa.pl
## USAGE: sbatch slurm_sucker_runbwa.sh

#SBATCH --job-name=runbwa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --account=ysctrout
#SBATCH --time=1-0:00:00
#SBATCH --mem=5G
#SBATCH -o /project/ysctrout/jjahner/suckers/slurms/std/stdout_runbwa
#SBATCH -e /project/ysctrout/jjahner/suckers/slurms/std/stderr_runbwa

module load arcc/1.0 gcc/12.2.0 bwa/0.7.17

cd /project/ysctrout/jjahner/suckers/

## new fastqs
## perl runbwa.pl fastqs/*.fastq

## reference individuals
perl runbwa.pl ref_forJosh/*.fastq
