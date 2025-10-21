# allelic_balance.R by SPJ 102025
# purpose: analyze AB on a per-site and per-ind basis as a final check before running
# sequoia on the full dataset (n = 1184)

# i decided i am going to start in sequoia with:
# path: /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem/vcfs/q40/mindep8/
# file: rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.recode.vcf

# all this comes from a dataset that has been...
# cleaned using illumina, phix, and ecoli databases, and has been cleaned using
# FASTP (remove bases with q<20, reads shorter than 50bp, poly g tails, and any 
# remaining sequence adapters). This dataset was then aligned to the walleye reference
# using BWA MEM, variant calling was run, vcf was reheadered, and the data were 
# filtered using the slurm_1 and slurm_3 filtering scripts for the following params.

# rawfiltered - see bcftools mpileup in slurm_sauger_variants.sh
# contam - contaminant filtering mentioned above
# fastp - see above
# bial - keep only biallelic sites
# no indels - keep sites that are not insertions or deletions
# q40 - keep sites with site quality > 40
# mindep8 - minimum mean read depth 4 or above (removes low quality sites, low confidence calls)
# maxdep75 - maximum mean read depth 75 or lower (removes paralogs)
# maf30 - minor allele frequency 30% or higher
# miss95 - each snp needs to have 5% or less missing data
# and thinned so that there are 100KB between each snp

# the output from the slurm scripts allelicBalance-Site.sh and allelicBalance-Samples.sh

library(dplyr)
setwd("/Users/samjohnson/Desktop/allelic_balance/")

# 943 snps
site_ab <- read.table(file = "meanPerSite_allelicBalance_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K_label.txt",
                      header = FALSE, sep = "\t") 

# setup/maintenence/good practice
colnames(site_ab) <- c("Scaffold", "Position", "AllelicBalance", "IgnoreLabel")
site_ab$AllelicBalance <-  as.numeric(site_ab$AllelicBalance)
head(site_ab)

# plot hist
hist(site_ab$AllelicBalance, breaks = 50, main = "AB Per Site at Heterozygous Sites (n snps = 943)",
     xlab = "AB for Heterozyotes")

# per-sample AB calculations had to be changed on the fly here. see second attempt
# in allelicBalance-Samples.sh
ind_ab <- read.table(file = "allelicBalance_rectangular_label.txt",
                     header = FALSE, sep = "\t")

# setup/maintenence/good practice
colnames(ind_ab) <- c("Scaffold", "Position", "IgnoreLabel", "Sample_ID", "AllelicBalance")
ind_ab$AllelicBalance <- as.numeric(ind_ab$AllelicBalance)
table(is.na(ind_ab$AllelicBalance)) # none!? i suppose the slurm filters out these options already with that if else statement

# create a new summary table, group by individual, and generate summary stats
summary_ind_ab <- ind_ab %>%
  group_by(Sample_ID) %>%
  summarize(
    n_sites = n(),               # number of heterozygous sites used
    mean_AB = mean(AllelicBalance, na.rm = TRUE),
    median_AB = median(AllelicBalance, na.rm = TRUE),
    sd_AB = sd(AllelicBalance, na.rm = TRUE)
  )

# convert back to df from tibble
summary_ind_ab <- data.frame(summary_ind_ab)

# plot hist
hist(summary_ind_ab$mean_AB, breaks = 50, main = "AB Per Sample at Heterozygous Sites (n inds = 1182)",
     xlab = "AB at Heterozygous Sites")
dim(summary_ind_ab)




