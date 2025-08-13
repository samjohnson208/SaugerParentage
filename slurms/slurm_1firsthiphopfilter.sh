#!/bin/sh

## slurm_1firsthiphopfilter.sh by SPJ 050525
## PURPOSE: filter the raw svit mem vcf to include biallelic sites, no indels, and reads with quality > 20
## USAGE: sbatch slurm_1firsthiphopfilter.sh

## note: also used 061525 to filter yellow perch vcfs for pca
## note: also used 062625 to filter yellow perch mem t2 vcfs for randomforest
## note: also used 070125 to filter walleye from the rehead_variants_rawfiltered_svit_mem_031325.vcf
## note: also used 070125 to filter yellow perch from the rehead_variants_rawfiltered_pflav_mem_031325.vcf

#SBATCH --job-name=firstfilter
#SBATCH --account=ysctrout
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_filters
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_filters
#SBATCH --mail-type=END

module load arcc/1.0 gcc/14.2.0 vcftools/0.1.17

# cd /project/ysctrout/hatchsauger/sam_sai_pflav_mem_t2/rf_2
cd /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem/vcfs

# filter raw vcf for biallelic sites
# vcftools --vcf rehead_variants_rawfiltered_pflav_mem_t2_062525.vcf --min-alleles 2 --max-alleles 2 --out variants_pflav_mem_t2_bial --recode  
#vcftools --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_080925.vcf --min-alleles 2 --max-alleles 2 --out rehead_variants_rawfiltered_svit_mem_contam_fastp_bial --recode  

# filter bial vcf for indels
# vcftools --vcf variants_pflav_mem_t2_bial.recode.vcf --remove-indels  --out variants_pflav_mem_t2_bial_noindels --recode
#vcftools --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial.recode.vcf --remove-indels --out rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels --recode

# now biallelic sites, no indels ----->  quality > 20
# vcftools --vcf variants_pflav_mem_t2_bial_noindels.recode.vcf --minQ 20 --out variants_pflav_mem_t2_bial_noindels_q20 --recode
#vcftools --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels.recode.vcf --minQ 20 --out rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20 --recode
#vcftools --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels.recode.vcf --minQ 30 --out rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q30 --recode
vcftools --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels.recode.vcf --minQ 40 --out rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40 --recode

