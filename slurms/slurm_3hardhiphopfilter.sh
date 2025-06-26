#!/bin/sh

## slurm_3hardhiphopfilter.sh by SPJ 050525 | modified 052725 for use to troubleshoot sequoia
## PURPOSE: filter the output from slurm_1firsthiphopfilter.sh according to a hard set of criteria (min/max depth, maf, and miss)
## USAGE: sbatch slurm_3hardhiphopfilter.sh

## note: also used 061525 to filter yellow perch vcfs for pca
## note: also used 062625 to filter yellow perch mem t2 vcfs for randomforest

#SBATCH --job-name=hardfilter
#SBATCH --account=ysctrout
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_pcahardfilter
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_pcahardfilter
#SBATCH --mail-type=END

module load arcc/1.0 gcc/14.2.0 vcftools/0.1.17

cd /project/ysctrout/hatchsauger/sam_sai_pflav_mem_t2

# filter first output file (variants_bial_noindels_q20.recode.vcf) by min and max mean depth per site across all samples.
vcftools --vcf variants_pflav_mem_t2_bial_noindels_q20.recode.vcf  --min-meanDP 6 --max-meanDP 75 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep8_maxdep75 --recode

# filter the vcf that now includes depth by maf and missing data per site.
vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep8_maxdep75.recode.vcf  --maf 0.01 --max-missing 0.9 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep6_maxdep75_maf1_miss90 --recode

# create genotype matrix
vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep8_maxdep75_maf1_miss95.recode.vcf --012 --out hard_variants_pflav_bial_mem_t2_noindels_q20_mindep8_maxdep75_maf1_miss95