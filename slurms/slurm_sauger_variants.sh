#!/bin/bash

## slurm_sauger_variants.sh by JPJ and SPJ 012325
## PURPOSE: to call variants
## USAGE: sbatch slurm_sauger_variants.sh

#SBATCH --job-name=variants_svit_mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=3-0:00:00
#SBATCH --mem=1G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_variants_svit_mem
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_variants_svit_mem

module load arcc/1.0 gcc/14.2.0 samtools/1.20 bcftools/1.20

cd /project/ysctrout/hatchsauger/sam_sai_pflav_mem_t2/bams/sorted.bams

# bcftools mpileup -a DP,AD,INFO/AD -C 50 -d 250 -f /project/ysctrout/reference_genomes/Sander_vitreus/walleye_genome.fna -q 30 -Q 20 -I -b bam_list.txt | bcftools call -v -m -f GQ -O z -o variants_rawfiltered_svit_mem_031325.vcf
bcftools mpileup -a DP,AD,INFO/AD -C 50 -d 250 -f /project/ysctrout/hatchsauger/Perca_flavescens/yellowperch.fna -q 30 -Q 20 -I -b bam_list.txt | bcftools call -v -m -f GQ -O z -o variants_rawfiltered_pflav_mem_t2_062525.vcf

