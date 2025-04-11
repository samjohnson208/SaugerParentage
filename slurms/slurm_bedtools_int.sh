#!/bin/bash

## slurm_bedtools_int.sh by SPJ 041125
## PURPOSE: to check the overlap in loci among two files
#3 USAGE: sbatch slurm_bedtools_int.sh

#SBATCH --job-name=bedtools_intersect
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --mem=1G
#SBATCH --time=06:00:00
#SBATCH --mem=4G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_bedtools_int
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_bedtools_int

module load arcc/1.0 gcc/14.2.0 samtools/1.20 bcftools/1.20

export PATH=/project/ysctrout/software/bedtools2/bin:$PATH

# define input files
FILE_A=/project/ysctrout/hatchsauger/sam_sai_svit_mem/rehead_variants_rawfiltered_svit_mem_031325.vcf
FILE_B=/project/ysctrout/hatchsauger/sam_sai_svit/rehead_variants_rawfiltered_svit_020625.vcf

# define output file
OUTPUT=/project/ysctrout/hatchsauger/bedtools_out/intersect_output_1.txt

# run intersect with count option
bedtools intersect -a "$FILE_A" -b "$FILE_B" -c > "$OUTPUT"