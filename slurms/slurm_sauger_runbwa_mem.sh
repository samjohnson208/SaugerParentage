#!/bin/bash

## slurm_sauger_runbwa.sh by JPJ and SPJ 112724
## PURPOSE: to use runbwa.pl
## USAGE: sbatch slurm_sauger_runbwa_mem.sh

#SBATCH --job-name=runbwa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --account=ysctrout
#SBATCH --time=3-0:00:00
#SBATCH --mem=5G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_runbwa_mem
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_runbwa_mem

module load arcc/1.0 gcc/14.2.0 bwa/0.7.17

cd /project/ysctrout/hatchsauger/

perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/runbwa_mem.pl /project/ysctrout/hatchsauger/1Saug/rawfastqs/*.fastq
