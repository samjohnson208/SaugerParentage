#!/bin/bash

## slurm_sucker_variants.sh by JPJ 7 vi 23
## PURPOSE: to call variants
## USAGE: sbatch slurm_sucker_variants.sh

#SBATCH --job-name=variants
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=0-12:00:00
#SBATCH --mem=1G
#SBATCH -o /project/ysctrout/jjahner/suckers/slurms/std/stdout_variants
#SBATCH -e /project/ysctrout/jjahner/suckers/slurms/std/stderr_variants

module load arcc/1.0 gcc/12.2.0 samtools/1.16.1 bcftools/1.16

cd /project/ysctrout/jjahner/suckers/sam_sai/

bcftools mpileup -a DP,AD,INFO/AD -C 50 -d 250 -f /project/ysctrout/reference_genomes/c_latipinnis/sucker_genome.fna -q 30 -Q 20 -I -b bam_list.txt | bcftools call -v -m -f GQ -O z -o variants_rawfiltered_11jun23.vcf


