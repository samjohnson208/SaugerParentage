#!/bin/bash

# allelicBalance-Site.sh by MPR on 072525
# purpose: to investigate allele balances for heterozygotes on a per-site basis

# Usage: sbatch this_script.sh input.vcf.gz label
# usage SPJ: sbatch /project/ysctrout/hatchsauger/SaugerParentage/slurms/allelicBalance-Site.sh /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem/vcfs/q40/mindep8/rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.recode.vcf label
# (from a directory where the results will be stored, in this case: /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem/vcfs/q40/mindep8/allelic_balance)
# the label is there just for a placeholder. MPR had a label in this script to visualize differences among libraries
# i don't have a use for that, but i'm not going to mess with changing the arch. of this script.

#SBATCH --job-name=ABhets-site
#SBATCH -A ysctrout
#SBATCH -t 0-05:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10   
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sjohn208@uwyo.edu


module spider | head
module load arcc/1.0 gcc/14.2.0 bcftools/1.20
 

vcf_file="$1"
label="$2"

outfile="meanPerSite_allelicBalance_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K_${label}.txt"

bcftools query -f '%CHROM\t%POS[\t%SAMPLE=%GT:%AD]\n' "$vcf_file" | \
awk -v tag="$label" '{
    chrom = $1;
    pos = $2;
    sum_ratio = 0;
    count = 0;

    for (i = 3; i <= NF; i++) {
        split($i, a, "=");
        split(a[2], b, ":");
        if (b[1] == "0/1" && b[2] ~ /^[0-9]+,[0-9]+$/) {
            split(b[2], ad, ",");
            ref = ad[1] + 0;
            alt = ad[2] + 0;
            if (ref + alt > 0) {
                ratio = alt / (ref + alt);
            } else if (alt > 0 && ref == 0) {
                ratio = 1;
            } else {
                ratio = 0;
            }
            sum_ratio += ratio;
            count++;
        }
    }

    if (count > 0) {
        mean_ratio = sum_ratio / count;
        printf "%s\t%s\t%.4f\t%s\n", chrom, pos, mean_ratio, tag;
    }
}' > "$outfile"






























































































































