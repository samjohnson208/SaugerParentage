#!/bin/sh

## slurm_bedLD by SPJ 090225
## PURPOSE: run LD-based site pruning using PLINK, then retain original VCF info
## USAGE: sbatch slurm_bedLD.sh

#SBATCH --job-name=bedLD
#SBATCH --account=ysctrout
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_bedLD
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_bedLD
#SBATCH --mail-type=END

module load arcc/1.0 gcc/14.2.0 vcftools/0.1.17

cd /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem/vcfs/sitequal/plink_LD

# step 1: convert VCF to PLINK format, set SNP IDs to CHR:POS
/project/ysctrout/mrodri23/programs/plink \
  --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_GQ20_mindep4_maxdep75_maf30_miss90.recode.vcf \
  --make-bed \
  --out bed_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_GQ20_mindep4_maxdep75_maf30_miss90 \
  --double-id \
  --allow-extra-chr \
  --set-missing-var-ids @:#

# step 2: LD-based pruning
/project/ysctrout/mrodri23/programs/plink \
  --bfile bed_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_GQ20_mindep4_maxdep75_maf30_miss90 \
  --indep-pairwise 5 2 0.2 \
  --allow-extra-chr \
  --out pruned_bed_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_GQ20_mindep4_maxdep75_maf30_miss90

# step 3a: create positions file (CHR and POS) from pruned variants
awk -F':' '{print $1"\t"$2}' \
  pruned_bed_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_GQ20_mindep4_maxdep75_maf30_miss90.prune.in \
  > pruned_positions.txt

# step 3b: filter original VCF using positions, retain all INFO and genotypes
vcftools \
  --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_GQ20_mindep4_maxdep75_maf30_miss90.recode.vcf \
  --positions pruned_positions.txt \
  --recode \
  --recode-INFO-all \
  --out pruned_bed_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_GQ20_mindep4_maxdep75_maf30_miss90





