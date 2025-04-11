#!/bin/bash

## slurm_bedtools_int.sh by SPJ 041125
## PURPOSE: to check the overlap in loci among two files
## USAGE: sbatch slurm_bedtools_int.sh

#SBATCH --job-name=bedtools_intersect
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_bedtools_int
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_bedtools_int

module load arcc/1.0 gcc/14.2.0 samtools/1.20 bcftools/1.20

export PATH=/project/ysctrout/software/bedtools2/bin:$PATH

# define input files
FILE_A=/project/ysctrout/hatchsauger/bedtools/svit_mem.bed
FILE_B=/project/ysctrout/hatchsauger/bedtools/svit_alnsamse.bed

# define output file
OUTPUT=/project/ysctrout/hatchsauger/bedtools/intersect_output_1.bed

# run intersect with count option
bedtools intersect -a "$FILE_A" -b "$FILE_B" > "$OUTPUT"
bedtools intersect -a file_a.bed -b file_b.bed -wa -wb > "$OUTPUT"
