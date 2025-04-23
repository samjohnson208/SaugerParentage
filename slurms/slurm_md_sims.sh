#!/bin/bash

## created by SPJ 0423
## PURPOSE: check tpfp for simulations on parent-offspring duos
## USAGE: sbatch slurm_md_sims.sh

#SBATCH --job-name=mdsims
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_mdsims
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_mdsims
#SBATCH --mail-type=END

# cd to directory with all of the pars files for each nloci, wild_samp, and md
/project/ysctrout/hatchsauger/SaugerParentage/simulations

# load modules
module load arcc/1.0 gcc/14.2.0 r/4.4.0

# run the R script!
Rscript sauger_cross_check_pars_sj.R





