#!/bin/sh

## slurm_sauger_permd.sh by SPJ on 021025
## PURPOSE: to create a .csv that lists the percent of missing data per site in a .vcf
## USAGE: sbatch slurm_sauger_permd.sh

#SBATCH --job-name=permissingdata
#SBATCH --account=ysctrout
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --mail-type=END
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_permd
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_permd

module load arcc/1.0
module load vcftools/0.1.17

cd /project/ysctrout/hatchsauger/sam_sai_pflav

vcftools --vcf rehead_variants_rawfiltered_pflav_012325.vcf --missing-site --out per_md_per_site_pflav

cd /project/ysctrout/hatchsauger/sam_sai_svit

vcftools --vcf rehead_variants_rawfiltered_svit_020625.vcf --missing-site --out per_md_per_site_svit



