#!/bin/bash

## created by SPJ 022525 | Modified by SPJ 032625
## PURPOSE: run sequoia on R on the cluster (on a converted genotype matrix)
## USAGE: sbatch slurm_sequoia_first.sh

#SBATCH --job-name=sequoia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_sequoia_first
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_sequoia_first
#SBATCH --mail-type=END

# cd to directory with your converted genotype matrix
cd /project/ysctrout/hatchsauger/SaugerParentage/r_scripts

# load modules
module load arcc/1.0 gcc/14.2.0 r/4.4.0

# run the R script!
Rscript Sequoia_Full.R





