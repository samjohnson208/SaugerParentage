#!/bin/sh

## slurm_3hardhiphopfilter.sh by SPJ 050525
## PURPOSE: filter the output from slurm_1firsthiphopfilter.sh according to a hard set of criteria (min/max depth, maf, and miss)
## USAGE: sbatch slurm_3hardhiphopfilter.sh

#SBATCH --job-name=hardfilter
#SBATCH --account=ysctrout
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_hiphophardfilter
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_hiphophardfilter
#SBATCH --mail-type=END

module load arcc/1.0 gcc/14.2.0 vcftools/0.1.17

cd /project/ysctrout/hatchsauger/sam_sai_svit_mem/HipHop_Inp

# filter first output file (variants_bial_noindels_q20.recode.vcf) by min and max mean depth per site across all samples.
vcftools --vcf variants_bial_noindels_q20.recode.vcf  --min-meanDP 8 --max-meanDP 75 --out hard_variants_bial_noindels_q20_mindep8_maxdep75 --recode

# filter the vcf that now includes depth by maf and missing data per site.
vcftools --vcf hard_variants_bial_noindels_q20_mindep8_maxdep75.recode.vcf  --maf 0.05 --max-missing 0.95 --out hard_variants_bial_noindels_q20_mindep8_maxdep75_maf5_miss95 --recode

# create genotype matrix
vcftools --vcf hard_variants_bial_noindels_q20_mindep8_maxdep75_maf5_miss95.recode.vcf --012 --out hard_variants_bial_noindels_q20_mindep8_maxdep75_maf5_miss95