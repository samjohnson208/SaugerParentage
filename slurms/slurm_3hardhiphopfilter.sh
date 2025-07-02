#!/bin/sh

## slurm_3hardhiphopfilter.sh by SPJ 050525 | modified 052725 for use to troubleshoot sequoia
## PURPOSE: filter the output from slurm_1firsthiphopfilter.sh according to a hard set of criteria (min/max depth, maf, and miss)
## USAGE: sbatch slurm_3hardhiphopfilter.sh

## note: also used 061525 to filter yellow perch vcfs for pca
## note: also used 062625 to filter yellow perch mem t2 vcfs for randomforest
## note: increased min depth to 15 on sequoia input (svit_mem maf 30 miss 95) on 062725
## note: also used 070125 to filter walleye from the rehead_variants_rawfiltered_svit_mem_031325.vcf
## note: also used 070125 to filter yellow perch from the rehead_variants_rawfiltered_pflav_mem_031325.vcf

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

cd /project/ysctrout/hatchsauger/sam_sai_pflav_mem_t2/rf_2

# filter first output file (variants_bial_noindels_q20.recode.vcf) by min and max mean depth per site across all samples.
# vcftools --vcf variants_pflav_mem_t2_bial_noindels_q20.recode.vcf  --min-meanDP 4 --max-meanDP 75 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75 --recode
# vcftools --vcf variants_pflav_mem_t2_bial_noindels_q20.recode.vcf  --min-meanDP 3 --max-meanDP 75 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75 --recode
# vcftools --vcf variants_bial_noindels_q20.recode.vcf --min-meanDP 15 --max-meanDP 75 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep15_maxdep75 --recode
# vcftools --vcf variants_pflav_mem_bial_noindels_q20.recode.vcf --min-meanDP 4 --max-meanDP 75 --out hard_variants_pflav_mem_bial_noindels_q20_mindep4_maxdep75 --recode
# vcftools --vcf variants_pflav_mem_t2_bial_noindels_q20.recode.vcf --min-meanDP 4 --max-meanDP 75 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75 --recode
# vcftools --vcf variants_pflav_mem_t2_bial_noindels_q20.recode.vcf --min-meanDP 3 --max-meanDP 75 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75 --recode

# filter the vcf that now includes depth by maf and missing data per site.
# vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75.recode.vcf  --maf 0.01 --max-missing 0.8 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75_maf1_miss80 --recode
# vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75.recode.vcf  --maf 0.01 --max-missing 0.8 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss80 --recode
# vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep15_maxdep75.recode.vcf --maf 0.30 --max-missing 0.95 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep15_maxdep75_maf30_miss95 --recode
# vcftools --vcf hard_variants_pflav_mem_bial_noindels_q20_mindep4_maxdep75.recode.vcf --maf 0.30 --max-missing 0.95 --out hard_variants_pflav_mem_bial_noindels_q20_mindep4_maxdep75_maf30_miss95 --recode
vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75.recode.vcf --maf 0.01 --max-missing 0.80 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75_maf1_miss80 --recode
vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75.recode.vcf --maf 0.01 --max-missing 0.80 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss80 --recode

# thin!
vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75_maf1_miss80.recode.vcf --thin 1000 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75_maf1_miss80_thin1K --recode
vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75_maf1_miss80.recode.vcf --thin 5000 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75_maf1_miss80_thin5K --recode
vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75_maf1_miss80.recode.vcf --thin 10000 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75_maf1_miss80_thin10K --recode

vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss80.recode.vcf --thin 1000 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss80_thin1K --recode
vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss80.recode.vcf --thin 5000 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss80_thin5K --recode
vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss80.recode.vcf --thin 10000 --out hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss80_thin10K --recode

# create genotype matrix
# vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep4_maxdep75_maf1_miss80.recode.vcf --012 --out hard_variants_pflav_bial_mem_t2_noindels_q20_mindep4_maxdep75_maf1_miss80
# vcftools --vcf hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss80.recode.vcf --012 --out hard_variants_pflav_bial_mem_t2_noindels_q20_mindep3_maxdep75_maf1_miss80



