#!/bin/sh

## slurm_1firsthiphopfilter.sh by SPJ 050525
## PURPOSE: filter the raw svit mem vcf to include biallelic sites, no indels, and reads with quality > 20
## USAGE: sbatch slurm_1firsthiphopfilter.sh

#SBATCH --job-name=firstfilter
#SBATCH --account=ysctrout
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_hiphop
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_hiphop
#SBATCH --mail-type=END

module load arcc/1.0 gcc/14.2.0 vcftools/0.1.17

cd /project/ysctrout/hatchsauger/sam_sai_svit_mem/HipHop_Inp

# filter raw vcf for biallelic sites
vcftools --vcf rehead_variants_rawfiltered_svit_mem_031325.vcf --min-alleles 2 --max-alleles 2 --out variants_bial --recode 

# filter bial vcf for indels
vcftools --vcf variants_bial.recode.vcf --remove-indels  --out variants_bial_noindels --recode

# now biallelic sites, no indels ----->  quality > 20
vcftools --vcf variants_bial_noindels.recode.vcf --minQ 20 --out variants_bial_noindels_q20 --recode




