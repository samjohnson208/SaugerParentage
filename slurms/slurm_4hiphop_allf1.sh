#!/bin/bash

## created by SPJ 050925
## PURPOSE: run hiphop on R on the cluster (on a converted genotype matrix)
## USAGE: sbatch slurm_hiphop_allf1

#SBATCH --job-name=hiphop_allf1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_hiphop_allf1
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_hiphop_allf1
#SBATCH --mail-type=END

# cd to directory with your converted genotype matrix and individuals file
cd /project/ysctrout/hatchsauger/sam_sai_svit_mem/HipHop_Inp/HipHop_AllF1

# load modules
module load arcc/1.0 gcc/14.2.0 r/4.4.0

# run the R script!
Rscript sauger_hiphop_allf1.R



