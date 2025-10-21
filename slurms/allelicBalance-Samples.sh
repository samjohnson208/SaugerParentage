#!/bin/bash

# allelicBalance-Site.sh by MPR on 072825
# purpose: to investigate allele balances for heterozygotes on a per-sample basis

# Usage: sbatch this_script.sh input.vcf.gz label
# usage SPJ: sbatch /project/ysctrout/hatchsauger/SaugerParentage/slurms/allelicBalance-Samples.sh /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem/vcfs/q40/mindep8/rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.recode.vcf label
# (from a directory where the results will be stored, in this case: /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem/vcfs/q40/mindep8/allelic_balance)
# the label is there just for a placeholder. MPR had a label in this script to visualize differences among libraries
# i don't have a use for that, but i'm not going to mess with changing the arch. of this script.


#SBATCH --job-name=ABhets-sample
#SBATCH -A ysctrout
#SBATCH -t 0-05:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10   
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sjohn208@uwyo.edu


module spider | head
#module load arcc/1.0 gcc/14.2.0 bcftools/1.20

#vcf_file="$1"
#label="$2"

#outfile="meanPerSample_allelicBalance_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K_label.txt"


#bcftools query -f '%CHROM\t%POS[\t%SAMPLE=%GT:%AD]\n' "$vcf_file" | awk -v lib="$label" '{
    #out = $1 "\t" $2 "\t" lib;
    #count = 0;

    #for (i=3; i<=NF; i++) {
        #split($i, a, "=");          # a[1] = sample name, a[2] = GT:AD
        #split(a[2], b, ":");        # b[1] = GT, b[2] = AD
        #if (b[1] == "0/1" && b[2] ~ /^[0-9]+,[0-9]+$/) {
            #split(b[2], ad, ",");   # ad[1] = ref, ad[2] = alt
            #ref = ad[1] + 0;
            #alt = ad[2] + 0;
            #if (ref + alt > 0) {
                #ratio = alt / (ref + alt);
            #} else if (alt > 0 && ref == 0) {
                #ratio = 1;
            #} else {
                #ratio = 0;
            #}
            #out = out "\t" a[1] "=" sprintf("%.4f", ratio);
            #count++;
        #}
    #}

    #if (count > 0) print out;
#}' > "$outfile"


# TROUBLESHOOTING.
# the last output isn't rectangular because not all samples are genotyped at every locus.
# i.e., each row is a different length.

module load arcc/1.0 gcc/14.2.0 bcftools/1.20

vcf_file="$1"
label="$2"

outfile="allelicBalance_rectangular_${label}.txt"

# Output format: CHROM  POS  LABEL  SAMPLE  RATIO
bcftools query -f '%CHROM\t%POS[\t%SAMPLE=%GT:%AD]\n' "$vcf_file" | \
awk -v lib="$label" '{
    chrom = $1;
    pos = $2;

    for (i=3; i<=NF; i++) {
        split($i, a, "=");          # a[1] = sample name, a[2] = GT:AD
        split(a[2], b, ":");        # b[1] = GT, b[2] = AD
        if (b[1] == "0/1" && b[2] ~ /^[0-9]+,[0-9]+$/) {
            split(b[2], ad, ",");   # ad[1] = ref, ad[2] = alt
            ref = ad[1] + 0;
            alt = ad[2] + 0;
            total = ref + alt;
            if (total > 0) {
                ratio = alt / total;
                printf "%s\t%s\t%s\t%s\t%.4f\n", chrom, pos, lib, a[1], ratio;
            }
        }
    }
}' > "$outfile"

# usage: sbatch /project/ysctrout/hatchsauger/SaugerParentage/slurms/allelicBalance-Samples.sh /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem/vcfs/q40/mindep8/rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.recode.vcf label

















































































































