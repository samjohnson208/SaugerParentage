#!/bin/sh

## slurm_rsq_linkage.sh by SPJ 082725
## PURPOSE: obtain correlation coefficients for pairwise combinations of snps in a vcf
## USAGE: sbatch slurm_rsq_linkage.sh

#SBATCH --job-name=LD_info
#SBATCH --account=ysctrout
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_rsq_linkage
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_rsq_linkage
#SBATCH --mail-type=END

module load arcc/1.0 gcc/14.2.0 vcftools/0.1.17

cd /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem/vcfs/q40/mindep8

vcftools --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95.recode.vcf --geno-r2  --ld-window-bp 100000 --out nothin_window100K
vcftools --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95.recode.vcf --geno-r2  --ld-window-bp 500000 --out nothin_window500K
vcftools --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin1M.recode.vcf --geno-r2 --ld-window-bp 1000000 --out thin1M_window1M
vcftools --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin2M.recode.vcf --geno-r2 --ld-window-bp 2000000 --out thin2M_window2M
vcftools --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin3M.recode.vcf --geno-r2 --ld-window-bp 3000000 --out thin3M_window3M
vcftools --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin4M.recode.vcf --geno-r2 --ld-window-bp 4000000 --out thin4M_window4M
vcftools --vcf rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin5M.recode.vcf --geno-r2 --ld-window-bp 5000000 --out thin5M_window5M


