## Sequoia_Contam_GQ_LD.R by SPJ 090225
## PURPOSE: Alright, here's the deal. We first wanted to see if we weren't getting
# many assignments due to genotype quality issues. We pushed minQ up to 4 and min
# mean depth up to 8, maf30 and miss95, but this still didn't get us to per indiv.
# site quality, nor did we have an idea of linkage with any statistical power. we
# used the distance based thinning in vcftools (--thin) but as you can see at the
# end of Sequoia_Contamfilt_mindep8.R, there are instances where there are sites
# that are close to one another but not linked, and ones far away from one another
# that are linked. I decided to use the --minGQ option and the LD filter implemented
# in plink to try to get some good sites that are not linked. what's strange about
# that r2 filter is that you set a sliding window based off of the number of snps, 
# rather than the number of bp, and the step for that sliding window is also based
# off of the number of snps rather than number of bp, so we really don't know if
# we're eliminating snps that are inherited together but far away from one another 
# physically. (high r2, high dist)

# that said, this dataset has been filtered as follows:

# rawfiltered: see bcftools mpileup in slurm_sauger_variants.sh
# contam: contaminant filtering (illumina, phix, ecoli databases)
# fastp: (remove bases with q<20, reads shorter than 50bp, poly g tails, and any 
  # remaining sequence adapters)
# bial: biallelic sites only
# noindels: no insertions/deletions
# q20: site quality (phred score based metric to detect if it is a TRUE variant site)
# GQ20: per individual, per site, genotype quality (phred score based)
# mindep4: site min mean read depth > 4 (removes low quality sites, low confidence calls)
# maxdep75: site max mean read depth < 75 (removes paralogs)
# maf30 - minor allele frequency 30% or higher
# miss90 - each snp needs to have 10% or less missing data

# additional notes: that GQ filter turns genotypes less than 20 into NA, so when
# you add that missing data filter, you just SLAUGHTER all of your sites if there
# aren't quality genotypes. that's why we had to go with miss 90. we can try to
# not do genotype qual but do --geno-r2 on this one:
# /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem/vcfs/q40/mindep8/
# rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95.recode.vcf

## libraries
################################################################################

# install.packages("sequoia")
# 65
library(sequoia)

# install.packages("dplyr")
# 65
library(dplyr)

################################################################################

setwd("/Users/samjohnson/Documents/Sauger_042225/GeneticData/Sequoia/Sequoia_Inp/contam_fastp_svit_mem/firstfilt_GQ_harddepth_hardmafmiss_LDpruned")

################################################################################
#### RUNNING ON FURTHER THINNED (LD FILTERED) FILES ####
################################################################################

ld_pruned_mat <- read.table(file = "pruned_bed_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_GQ20_mindep4_maxdep75_maf30_miss90.012_conv",
                            header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

# remove first column. ascending numbers 1 - 1184
ld_pruned_mat <- ld_pruned_mat[, -1]

# as matrix
ld_pruned_mat <- as.matrix(ld_pruned_mat)

# yes, i am aware that these are all the same. trying to be thorough. would hate 
# to mess up a matrix because i was lazy on a step like this.
ind_ld_pruned_mat <- read.table(file = "pruned_bed_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_GQ20_mindep4_maxdep75_maf30_miss90.012.indv", header = FALSE)

# rename first column to something meaningful
ind_ld_pruned_mat <- ind_ld_pruned_mat  %>% 
  rename(sample = V1)

# establish sample identities in the geno_mat
rownames(ld_pruned_mat) <- ind_ld_pruned_mat$sample

# read in scaffolds and positions
pos_ld_pruned_mat <- read.table(file = "pruned_bed_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_GQ20_mindep4_maxdep75_maf30_miss90.012.pos", header = FALSE)
pos_ld_pruned_mat$position <- paste(pos_ld_pruned_mat$V1, pos_ld_pruned_mat$V2, sep = "_")

# set positions as column names of geno_mat
colnames(ld_pruned_mat) <- pos_ld_pruned_mat$position

# read in parent offspring (F0 and Test F1 sample id's)
testsamp <- read.csv(file = "posampleids.csv", header = TRUE) # 210 samples, 114 Parents, 96 Test F1's

# read in lh data
LH_Data <- read.csv(file = "testindivs_LH.csv", header = TRUE)
# alright, what's going on here with the missing individual?

missingind <- testsamp$sample[!testsamp$sample %in% LH_Data$ID]
# first thing is to check the extractions, readme, and hiphop scripts.
# i understand 6757 wasn't sequenced, but what's up with 6436. it's on the plate map!
# 6436 wasn't spawned. guess i didn't catch that when i was prepping for sequencing.
# test group is for sure 208 individuals. 113 F0 parents, 95 test F1's.

# change the F's to 1 and M's to 2
LH_Data$Sex[LH_Data$Sex == "M"] <- 2
LH_Data$Sex[LH_Data$Sex == "F"] <- 1
LH_Data$Sex <- as.numeric(LH_Data$Sex)
str(LH_Data)
table(LH_Data$Sex)

# filter the genotype matrices so that it only includes the ids from LH_Data$ID
ld_pruned_mat <- ld_pruned_mat[rownames(ld_pruned_mat) %in% LH_Data$ID, , drop = FALSE]
dim(ld_pruned_mat) # 208 indivs, 259 loci

# check for all heterozygous sites since those likely will not be informative.
all_ones_ld_pruned_mat <- apply(ld_pruned_mat, 2, function(col) all(col == 1))
ld_pruned_mat <- ld_pruned_mat[, !all_ones_ld_pruned_mat]
dim(ld_pruned_mat) # 208 indivs, 258 loci, lost one all-heterozygote snp

# check genotype matrix for samples/loci to be excluded
check_ld_pruned_mat <- CheckGeno(ld_pruned_mat, quiet = FALSE, Plot = TRUE, 
                                 Return = "GenoM", Strict = TRUE, 
                                 DumPrefix = c("F0", "M0"))
# ✖ There are 1 individuals scored for <20% of SNPs,it is advised to treat their assignments with caution
# ℹ There are  208  individuals and  258  SNPs.

# who is that one indiv?
na_per_row <- rowSums(check_ld_pruned_mat == -9)
max_na <- max(na_per_row)
rows_most_na <- which(na_per_row == max_na)

# SAR_15_6437 is genotyped for less than 20% of the snps. remove it. 
check_ld_pruned_mat <- check_ld_pruned_mat[rownames(check_ld_pruned_mat) != "SAR_15_6437", ]
dim(check_ld_pruned_mat) # 207 indivs, 258 snps.

# check again
check_ld_pruned_mat <- CheckGeno(check_ld_pruned_mat, quiet = FALSE, Plot = TRUE, 
                                 Return = "GenoM", Strict = TRUE, 
                                 DumPrefix = c("F0", "M0"))
# ✔ Genotype matrix looks OK! There are  207  individuals and  258  SNPs.

# are there any sites with no missing data?
completesites_ld_pruned_mat <- apply(check_ld_pruned_mat, 2, function(col) all(col != -9))
table(completesites_ld_pruned_mat) # OUCH. ONLY TWO WITH NO MD.

erm1 <- ErrToM(Err = 0.01)
erm25 <- ErrToM(Err = 0.025)
erm5 <- ErrToM(Err = 0.05)

gmr_ld_pruned_mat_erm5 <- GetMaybeRel(GenoM = check_ld_pruned_mat,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_ld_pruned_mat))

length(unique(gmr_ld_pruned_mat_erm5[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_ld_pruned_mat_erm5[["MaybeTrio"]]$id)]))

# erm1: 123, 68, 105. 258 snps. 82 test f1's in the trios. 
# erm25: 121, 76, 149, 258 snps. 83 test f1's in the trios.
# erm5: 112, 144, 207, 258 snps. 82 test f1's in the trios. 

##### ----- Alright, why don't we then try to validate each of these ----- #####

trios_gmr_ld_pruned_erm5 <- gmr_ld_pruned_mat_erm5[["MaybeTrio"]]
head(trios_gmr_ld_pruned_erm5)

# read in lookup table
cross_lookup <- read.csv(file = "true_cross_lookup.csv", header = TRUE)
cross_lookup <- cross_lookup[,-1]

# this is a ridiculously effective set of piped functions here:
trios_gmr_ld_pruned_erm5 <- trios_gmr_ld_pruned_erm5 %>%
  # create a new column called pair, where the cross is structured the same way
  # as those in the lookup table. again, we're using pmin and pmax. sorting the
  # pair entries this way ensures that the cross A_B will be treated the same as
  # the cross B_A.
  mutate(pair = paste(pmin(parent1, parent2), pmax(parent1, parent2), sep = "_")) %>%
  
  left_join(cross_lookup %>% # join the top3 dataframe to the lookup table
              select(pair) %>% # but ONLY the column cross_lookup$pair
              distinct() %>% # keeps only unique entries of the known pairs, avoids any repeats
              mutate(valid_cross = TRUE), # creates a new column called valid_cross, which is either TRUE
            # or NA, if the pair exists in the lookup table, or if it does not.
            by = "pair") %>% # actually completes the join, by the shared pair column.
  # this is the line in which the matching actually occurs. 
  mutate(valid_cross = ifelse(is.na(valid_cross), FALSE, valid_cross)) # recreates that valid_cross column but turns
# all of the NA's formed by the join into FALSE's

head(trios_gmr_ld_pruned_erm5)
table(trios_gmr_ld_pruned_erm5$valid_cross)


trios_gmr_ld_pruned_erm5_checked <- trios_gmr_ld_pruned_erm5 %>%
  select(id, parent1, parent2, pair, valid_cross, everything())

trios_gmr_ld_pruned_erm5_checked <- trios_gmr_ld_pruned_erm5_checked %>% 
  filter(!LLRparent1 %in% c(555, -555),
         !LLRparent2 %in% c(555, -555),
         !LLRpair    %in% c(555, -555))
dim(trios_gmr_ld_pruned_erm5_checked)
table(trios_gmr_ld_pruned_erm5_checked$valid_cross)

length(unique(trios_gmr_ld_pruned_erm1_checked$id[grepl("^SAR_15_67", trios_gmr_ld_pruned_erm1_checked$id)]))
length(unique(trios_gmr_ld_pruned_erm25_checked$id[grepl("^SAR_15_67", trios_gmr_ld_pruned_erm25_checked$id)]))
length(unique(trios_gmr_ld_pruned_erm5_checked$id[grepl("^SAR_15_67", trios_gmr_ld_pruned_erm5_checked$id)]))

# erm1: 82 unique test f1's assigned, 86.32%, 72 true, 10 false. 87.80% 
# erm25: 83 unique test f1's assigned, 87.38%, 73 true, 11 false. 87.95%
# erm5: 82 unique test f1's assigned, 86.32%, 71, 14 false. 86.58%

# alright, well at least it's good to know that this LD and GQ filtering didn't 
# really improve anything. now we can focus efforts on finishing the table.


