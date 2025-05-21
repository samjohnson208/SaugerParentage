#!/bin/sh

## slurm_2softhiphopfilter.sh by SPJ 050525
## PURPOSE: filter the output from slurm_1firsthiphopfilter.sh according to a soft set of criteria (no depth, just maf and miss)
## USAGE: sbatch slurm_2softhiphopfilter.sh

#SBATCH --job-name=softfilter
#SBATCH --account=ysctrout
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_hiphopsoftfilter
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_hiphopsoftfilter
#SBATCH --mail-type=END

module load arcc/1.0 gcc/14.2.0 vcftools/0.1.17

cd /project/ysctrout/hatchsauger/sam_sai_svit_mem/HipHop_Inp

# filter first output file (variants_bial_noindels_q20.recode.vcf) by maf and missing data per site. no depth.
vcftools --vcf variants_bial_noindels_q20.recode.vcf  --maf 0.01 --max-missing 0.9 --out soft_variants_bial_noindels_q20_maf1_miss90 --recode

# create genotype matrix
vcftools --vcf soft_variants_bial_noindels_q20_maf1_miss90.recode.vcf --012 --out soft_variants_bial_noindels_q20_maf1_miss90

