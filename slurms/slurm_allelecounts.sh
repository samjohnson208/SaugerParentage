#!/bin/sh

## slurm_allelecounts.sh by SPJ and MPR 101025
## PURPOSE: count the number of each allele for each site in a vcf
## USAGE: sbatch slurm_allelecounts.sh

#SBATCH --job-name=allelecounts
#SBATCH --account=ysctrout
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_allelecounts
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_allelecounts
#SBATCH --mail-type=END

module load arcc/1.0  gcc/14.2.0 bcftools/1.20

{echo -e “chrom\tpos\tn_ref\tn_het\tn_alt\tn_miss”; bcftools query -f ‘%CHROM\t%POS[\t%GT]\n’ rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin1M.recode.vcf \ | awk ‘{n_ref=0; n_het=0; n_alt=0; n_miss=0; for (i=3; i<=NF; i++) {if ($i==“0/0" || $i==“0|0”) n_ref++; else if ($i==“1/1" || $i==“1|1”) n_alt++; else if ($i==“0/1" || $i==“1/0” || $i==“0|1" || $i==“1|0”) n_het++; else if ($i==“./.” || $i==“.|.“) n_miss++;} print $1, $2, n_ref, n_het, n_alt, n_miss;}’} | column -t > site_counts_thin1M.txt
	
	
	