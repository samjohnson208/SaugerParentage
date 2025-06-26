#!/bin/bash


#SBATCH --job-name=FiltSiteDepth
#SBATCH --account=ysctrout
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module load arcc/1.0 gcc/14.2.0 vcftools/0.1.17

vcftools --gzvcf variants_pflav_mem_t2_bial_noindels_q20.recode.vcf --site-mean-depth --out variants_pflav_mem_t2_bial_noindels_q20

