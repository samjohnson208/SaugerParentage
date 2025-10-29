# Sequoia_ContamFilt_mindep8_md5_RADseqErr_ALLINDS.R by SPJ 102225
## PURPOSE: I am ready to try GetMaybeRel() on all three generations of fish.
# I have decided on a RADseq error matrix and a dataset, here they are:

# path: /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem/vcfs/q40/mindep8/
# file: rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.recode.vcf

# the data have been...
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

# work 102825 - think the first step in analyzing the pedigree is to take the LLRs
# for both duos and trios assigned here, and map it against the dist of those for
# the test group. first step is to get both genomats in here together. starting here:

## libraries
################################################################################

# install.packages("sequoia")
# 65
library(sequoia)

# install.packages("dplyr")
# 65
library(dplyr)

################################################################################

setwd("/Users/samjohnson/Documents/Sauger_102325/GeneticData/Sequoia/Sequoia_Inp/contam_fastp_svit_mem/firstfilt_hardfilt_thinned/mindep8_maf30/geno_mat")

################################################################################


##### ---- Preparing GenoM ---- #####
mat_thin100K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012_conv",
                           header = FALSE, sep = "\t", na.strings = c("NA", "-1"))
# remove first column. ascending numbers 1 - 1184
mat_thin100K <- mat_thin100K[, -1]
# turn to matrix
mat_thin100K <- as.matrix(mat_thin100K)
dim(mat_thin100K) # 943 snps

# add inds
ind100K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012.indv", header = FALSE)
ind100K <-  ind100K %>% 
  rename(sample = V1)
rownames(mat_thin100K) <- ind100K$sample

# read in scaffolds and positions
pos100K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012.pos", header = FALSE)
dim(pos100K)

# create full positions by combining the two columns
pos100K$position <-  paste(pos100K$V1, pos100K$V2, sep = "_")

# set positions as column names of GenoM
colnames(mat_thin100K) <- pos100K$position

##### ----- ---- #####

##### ---- Adding Additional Data, Filtering by Individuals ---- #### 

# --- LH FOR ALL --- #
# first i need to add LH data, so that I can filter the GenoM to these indivs.
LH_All <- read.csv(file = "LH_F0_F1Spawn_F1Juv_F2.csv", header = TRUE)
# change the F's to 1 and M's to 2, all others are 3's
LH_All$Sex[LH_All$Sex == "M"] <- 2
LH_All$Sex[LH_All$Sex == "F"] <- 1
LH_All$Sex <- as.numeric(LH_All$Sex)
str(LH_All)
dim(LH_All)

table(LH_All$Sex)
# Both 15 and 16 F0's
# 134 males, 129 females

table(LH_All$Sex)
# Both 15 and 16 F0's
# 134 males, 129 females

table(LH_All$BirthYear)
# 263 F0's (15 and 16), 6436 was not spawned, so it's not included here.
# 334 F1's (juveniles from fall 15 and spawn agg fish from spring 21)
# NO TEST F1's included
# 480 F2's (19, 20, 21, all over)
# = 1077
# note: will's samples excluded, possible hybrids excluded, WF fish excluded
# THIS IS READY NOW FOR GETMAYBEREL



# --- LH FOR TEST --- #
# read in parent offspring (F0 and Test F1 sample id's)
testsamp <- read.csv(file = "posampleids.csv", header = TRUE) # 210 samples, 114 Parents, 96 Test F1's

# read in lh data for test inds
LH_Test <- read.csv(file = "testindivs_LH.csv", header = TRUE)
# alright, what's going on here with the missing individual?

missingind <- testsamp$sample[!testsamp$sample %in% LH_Test$ID]
# first thing is to check the extractions, readme, and hiphop scripts.
# i understand 6757 wasn't sequenced, but what's up with 6436. it's on the plate map!
# 6436 wasn't spawned. guess i didn't catch that when i was prepping for sequencing.
# test group is for sure 208 individuals. 113 F0 parents, 95 test F1's.

# change the F's to 1 and M's to 2
LH_Test$Sex[LH_Test$Sex == "M"] <- 2
LH_Test$Sex[LH_Test$Sex == "F"] <- 1
LH_Test$Sex <- as.numeric(LH_Test$Sex)
str(LH_Test)
table(LH_Test$Sex) # 53 female F0, 60 male F0, 95 unknown sex test F1
LH_Test$BirthYear[1:95] <- 1
LH_Test$BirthYear[96:nrow(LH_Test)] <- 0
# THIS IS READY NOW FOR GETMAYBEREL



# --- Filtering the mat_thin100K GenoM to include each of the groups of inds --- #

# filter GenoM so it only includes ids from LH_All$ID
dim(mat_thin100K) # 1184, 943
mat_thin100K_all <- mat_thin100K[rownames(mat_thin100K) %in% LH_All$ID, , drop = FALSE]
dim(mat_thin100K_all) # 1060 943
dim(LH_All)       # 1077 3

# where did those other 17 inds go?
missinginds <- LH_All$ID[!LH_All$ID %in% rownames(mat_thin100K)]
# SAR_21_5582 - missing fin clip
# SAR_21_5641 - missing fin clip
# SAR_21_5647 - missing fin clip
# SAR_21_5652 - missing fin clip
# SAR_21_5771 - missing fin clip
# SAR_21_5819 - missing fin clip
# SAR_19_5963 - missing fin clip
# SAR_19_5994 - missing fin clip
# SAR_21_6189 - lost fin clip
# SAR_21_6294 - missing fin clip
# SAR_21_6295 - missing fin clip
# SAR_21_6296 - missing fin clip
# SAR_21_6297 - missing fin clip
# SAR_21_6298 - missing fin clip
# SAR_21_6299 - missing vial
# SAR_21_6300 - missing vial
# SAR_21_6398 - missing fin clip

# okay, good to know i'm not an organizational disaster
# this happened because i pulled all of the sample id's from SAR_Data_092424 to 
# make this LH_All object. these are all of the tags/vials/clips i had rather than
# what actually got plated. 

# filter out missing inds
LH_All <- LH_All %>% 
    filter(ID %in% rownames(mat_thin100K_all))
dim(mat_thin100K_all) # 1060 943
dim(LH_All)       # 1060 3
# okay, that's remedied. great.

# filter GenoM so it only includes ids from LH_Test$ID
dim(mat_thin100K) # 1184, 943
mat_thin100K_test <- mat_thin100K[rownames(mat_thin100K) %in% LH_Test$ID, , drop = FALSE]
dim(mat_thin100K_test) # 208 943
dim(LH_Test)       # 208 3

# so now we're dealing with two pairs of objects: mat_thin100K_all with LH_All
                                           # and: mat_thin100K_test with LH_Test

dim(mat_thin100K_all)
dim(LH_All)

dim(mat_thin100K_test)
dim(LH_Test)

# nice.
##### ----- ---- #####

##### ---- HWE Filtering and Final GenoMat Maintenence ---- #####
# except no HWE filtering here.

# WORK ON FULL DATASET (mat_thin100K_all)
# check for all heterozygous sites since those likely will not be informative.
all_ones_100_all <- apply(mat_thin100K_all, 2, function(col) all(col == 1))
# remove all het. sites
mat_thin100K_all <- mat_thin100K_all[, !all_ones_100_all]

dim(mat_thin100K_all) # retained all

check_thin100K_all <- CheckGeno(mat_thin100K_all, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
# ✖ There are 2 individuals scored for <5% of SNPs, these WILL BE IGNORED
# ✖ In addition, there are 6 individuals scored for <20% of SNPs, it is advised to treat their assignments with caution
# ℹ After exclusion, There are  1058  individuals and  943  SNPs.

# should consider filtering for per individual md!!!
# let's just let it fly right now and see what happens. we'll tighten that up later.

# WORK ON TEST DATASET (mat_thin100K_test)
# check for all heterozygous sites since those likely will not be informative.
all_ones_100_test <- apply(mat_thin100K_test, 2, function(col) all(col == 1))
# remove all het. sites
mat_thin100K_test <- mat_thin100K_test[, !all_ones_100_test]
dim(mat_thin100K_test) # lost one

check_thin100K_test <- CheckGeno(mat_thin100K_test, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
# ✔ Genotype matrix looks OK! There are  208  individuals and  942  SNPs.
##### ----- ---- #####

##### ---- Missing Data per Individual ---- #####
# first runs say that there are six individuals genotyped for <20% of snps. 
# let's remedy that.

miss_per_ind <- read.table(file = "miss_per_indv_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.imiss", header = TRUE)
dim(miss_per_ind)
inds_to_keep <- miss_per_ind$INDV[miss_per_ind$F_MISS < 0.2]
inds_to_keep # individuals who are genotyped for 80% or samples or MORE
length(inds_to_keep) # 1154
dim(miss_per_ind) # 1184, means 30 inds are filtered out...

dim(check_thin100K_all) # 1058 samples from JUST the groups we're interested in (F0, F1Spawn, F1Juv, F2)
check_thin100K_all <- check_thin100K_all[rownames(check_thin100K_all) %in% inds_to_keep, ]
dim(check_thin100K_all) # 1030 samples remaining.

dim(check_thin100K_test) # 208 inds
check_thin100K_test <- check_thin100K_test[rownames(check_thin100K_test) %in% inds_to_keep, ]
dim(check_thin100K_test) # lost two inds... somehow...

dim(check_thin100K_all) # 1030 inds
dim(check_thin100K_test) # 206 inds


##### ---- ---- #####

# checking LLRs of the full dataset against those from the test. going into sequoia()
# and GetMaybeRel() here are what the matrices should look like. correspond to LH_All
# and LH_Test dataframes. These HAVE been filtered for md per individual.
dim(check_thin100K_all) # 1030 inds
dim(check_thin100K_test) # 206 inds

##### ---- sequoia() and GetMaybeRel() ---- #####

errM <- Err_RADseq(E0 = 0.075, E1 = 0.025, Return = 'matrix')
errM

seq <- sequoia(GenoM = check_thin100K, 
               LifeHistData = LH_All, 
               Module = "ped", 
               Err = errM,
               Complex = "full", 
               Herm = "no", 
               UseAge = "yes",
               StrictGenoCheck = TRUE, 
               DummyPrefix = c("F", "M"),
               Tfilter = -2, 
               Tassign = 1.0)
# ✔ assigned 112 dams and 105 sires to 1058 + 71 individuals (real + dummy)
# ✔ assigned 106 dams and 101 sires to 1030 + 66 individuals (real + dummy) 
# run two after removing individuals with 20% or more missing data. surprised
# that there were fewer assignments, though there are fewer inds i suppose...
# I hope that this also means that more of these assignments are accurate.


# overlapping gens? going to have to fix this. bet that ups the number of assignments.
# that's obviously the next step.

# 10/26/25, trying to mess with age stuff. let's see how we can assign parents with
# discrete generations on the most recent check_thin100K, which had inds with high md
# removed from it. dim(check_thin100K) = 1030 943
seq <- sequoia(GenoM = check_thin100K_10, 
               LifeHistData = LH_All, 
               Module = "ped",
               Err = errM,
               args.AP=list(Discrete = TRUE),
               Complex = "full", 
               Herm = "no", 
               UseAge = "yes", 
               StrictGenoCheck = TRUE, 
               DummyPrefix = c("F", "M"),
               Tfilter = -2, 
               Tassign = 1.0)
# with discrete generations:
# with module = "par": ✔ assigned 27 dams and 24 sires to 1030 individuals 
# with module = "ped": ✔ assigned 106 dams and 102 sires to 1030 + 71 individuals (real + dummy) 

# so now we need to start subsampling loci, see if that has an effect...
# subsample to 90%, 80%, 70%, so on...
samp_loci_90 <- sample(ncol(check_thin100K), 849)
check_thin100K_90 <- check_thin100K[, samp_loci_90]

samp_loci_80 <- sample(ncol(check_thin100K), 754)
check_thin100K_80 <- check_thin100K[, samp_loci_80]

samp_loci_70 <- sample(ncol(check_thin100K), 660)
check_thin100K_70 <- check_thin100K[, samp_loci_70]

samp_loci_60 <- sample(ncol(check_thin100K), 566)
check_thin100K_60 <- check_thin100K[, samp_loci_60]

samp_loci_50 <- sample(ncol(check_thin100K), 472)
check_thin100K_50 <- check_thin100K[, samp_loci_50]

samp_loci_40 <- sample(ncol(check_thin100K), 377)
check_thin100K_40 <- check_thin100K[, samp_loci_40]

samp_loci_30 <- sample(ncol(check_thin100K), 283)
check_thin100K_30 <- check_thin100K[, samp_loci_30]

samp_loci_20 <- sample(ncol(check_thin100K), 189)
check_thin100K_20 <- check_thin100K[, samp_loci_20]

samp_loci_10 <- sample(ncol(check_thin100K), 94)
check_thin100K_10 <- check_thin100K[, samp_loci_10]

dim(check_thin100K_90)
dim(check_thin100K_80)
dim(check_thin100K_70)
dim(check_thin100K_60)
dim(check_thin100K_50)
dim(check_thin100K_40)
dim(check_thin100K_30)
dim(check_thin100K_20)
dim(check_thin100K_10)

# alright, those are set. here are the results for each of those, structured this way:

# seq <- sequoia(GenoM = check_thin100K, 
#                LifeHistData = LH_All, 
#                Module = "ped",
#                Err = errM,
#                args.AP=list(Discrete = TRUE),
#                Complex = "full", 
#                Herm = "no", 
#                UseAge = "yes", 
#                StrictGenoCheck = TRUE, 
#                DummyPrefix = c("F", "M"),
#                Tfilter = -2, 
#                Tassign = 1.0)


# 90: ✔ assigned 106 dams and 103 sires to 1030 + 72 individuals (real + dummy)
# 80: ✔ assigned 106 dams and 97 sires to 1030 + 71 individuals (real + dummy) 
# 70: ✔ assigned 105 dams and 97 sires to 1030 + 69 individuals (real + dummy)
# 60: ✔ assigned 95 dams and 87 sires to 1030 + 60 individuals (real + dummy) 
# 50: ✔ assigned 101 dams and 93 sires to 1030 + 66 individuals (real + dummy) 
# 40: ✔ assigned 95 dams and 87 sires to 1030 + 61 individuals (real + dummy) 
# 30: ✔ assigned 88 dams and 86 sires to 1030 + 60 individuals (real + dummy) 
# 20: ✔ assigned 71 dams and 69 sires to 1030 + 51 individuals (real + dummy) 
# 10: ✔ assigned 44 dams and 44 sires to 1030 + 30 individuals (real + dummy) 

# consistently losing assignments with fewer loci. 


gmr <- GetMaybeRel(GenoM = check_thin100K,
                   Err = errM,
                   Module = "ped",
                   Complex = "full",
                   LifeHistData = LH_All,
                   quiet = TRUE,
                   Tfilter = -2,
                   Tassign = 1.0,
                   MaxPairs = 7 * nrow(check_thin1M))

# work 102725: I'm becoming very interested in using the LLRs to my advantage:

# reran today with 

# seq <- sequoia(GenoM = check_thin100K, 
#                LifeHistData = LH_All, 
#                Module = "ped", 
#                Err = errM,
#                Complex = "full", 
#                Herm = "no", 
#                UseAge = "yes",
#                CalcLLR = TRUE,
#                StrictGenoCheck = TRUE, 
#                DummyPrefix = c("F", "M"),
#                Tfilter = -2, 
#                Tassign = 1.0)

# ✔ assigned 106 dams and 101 sires to 1030 + 66 individuals (real + dummy) 

# why don't we then take these LLR's and see how they align with those from the
# test group. Should have an LLR for each. See how those compare for dummy, real
# inds for each gen/group. (e.g., 15 vs 16 parents, etc.), and just how many got
# assigned for each.

# work 102825: created the setup to run the pedigree reconstruction on the test
# set as well as the full set to compare LLRs.

# input object are as follows
dim(check_thin100K_all) # 1030 inds
dim(check_thin100K_test) # 206 inds
dim(LH_All) # 1060 inds (have not been filtered for md/ind... does this matter)
dim(LH_Test) # 208 inds (same here...)

seq_all <- sequoia(GenoM = check_thin100K_all,
                  LifeHistData = LH_All,
                  Module = "ped",
                  Err = errM,
                  Complex = "full",
                  Herm = "no",
                  UseAge = "yes",
                  CalcLLR = TRUE,
                  StrictGenoCheck = TRUE,
                  DummyPrefix = c("F", "M"),
                  Tfilter = -2,
                  Tassign = 1.0)
# ✔ assigned 106 dams and 101 sires to 1030 + 66 individuals (real + dummy)

seq_test <- sequoia(GenoM = check_thin100K_test,
                   LifeHistData = LH_Test,
                   Module = "ped",
                   Err = errM,
                   Complex = "simp",
                   Herm = "no",
                   UseAge = "yes",
                   CalcLLR = TRUE,
                   StrictGenoCheck = TRUE,
                   DummyPrefix = c("F", "M"),
                   Tfilter = -2,
                   Tassign = 1.0)
# Complex = full for both modules
# Module = "ped"
# ✔ assigned 94 dams and 95 sires to 206 + 3 individuals (real + dummy)
# Module = "par"
# ✔ assigned 90 dams and 90 sires to 206 individuals 

# Complex = "simp"
# Module = "par"
# ✔ assigned 91 dams and 90 sires to 206 individuals 

# not incredibly different, which is nice. goal here is to make comparisons

################################################################################
# after this, how different is the sequoia test object than the gmr object?

# investigate seq_test object (write code to investigate these objects)

# rerun gmr to get the same output from what we've gotten already

# run the valid_cross pipeline on that gmr object, generate vectors of the LLRs
  # for valid cross pairs and create a histogram of it.
################################################################################
# here's the code from Sequoia_ContamFilt_mindep8_md5_RADseqErr.R (slightly modified)
# e.g., included the Tfilter = -2 to match the sequoia runs above
# run gmr
gmr_test <- GetMaybeRel(GenoM = check_thin100K_test,                    # here (2)
                          Err = errM,                                   # error
                          # SeqList = outfull,
                          Module = "par",
                          # MaxMismatch = NA,
                          Complex = "simp",
                          LifeHistData = LH_Test,
                          quiet = FALSE,
                          Tassign = 1.0,
                          Tfilter = -2,
                          Herm = "no",
                          MaxPairs = 7*nrow(check_thin100K_test))      # here (1)
# ✔ Found 175 likely parent-offspring pairs, and 12, other non-assigned pairs of possible relatives
# ✔ Found 204 parent-parent-offspring trios


#unique test f1's in the trios:
length(unique(gmr_test[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_test[["MaybeTrio"]]$id)])) # here (2)
# 89

# how many unique focal indivs are in the trios?
length(unique(gmr_test[["MaybeTrio"]]$id))                             # here (1)
# 111

# see who it is
table(unique(gmr_test[["MaybeTrio"]]$id))                              # here (1)

trios_test <- gmr_test[["MaybeTrio"]]                                  # here (2)
head(trios_test)                                                       # here (1)

# read in lookup table
cross_lookup <- read.csv(file = "true_cross_lookup.csv", header = TRUE)
cross_lookup <- cross_lookup[,-1]

# this is a ridiculously effective set of piped functions here:
trios_test <- trios_test %>%                                            # here (2)
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
  mutate(valid_cross = ifelse(is.na(valid_cross), FALSE, valid_cross)) 
# recreates that valid_cross column but turns
# all of the NA's formed by the join into FALSE's

head(trios_test)                                                            # here (1)
dim(trios_test)                                                             # here (1)
table(trios_test$valid_cross)                                               # here (1)
# assign was 79/89 (88.76%), accuracy was 89/95 (93.68%)
# composite was 83.158 YES!!!

trios_test_checked <- trios_test %>%                                    # here (2)
  select(id, parent1, parent2, pair, valid_cross, everything())

trios_test_checked <- trios_test_checked %>%                            # here (2)
  filter(!LLRparent1 %in% c(555, -555),
         !LLRparent2 %in% c(555, -555),
         !LLRpair    %in% c(555, -555))
dim(trios_test_checked)                                                 # here (1)
# 89 checked trios
table(trios_test_checked$valid_cross)                                   # here (1)
# 79 valid trios
length(unique(trios_test_checked$id[grepl("^SAR_15_67", trios_test_checked$id)])) # here (2)
# 89 checked trios

# create a vector and histogram of the LLRs for each duo and trio
head(trios_test_checked)
table(is.na(trios_test_checked$LLRparent1)) # no na's
table(is.na(trios_test_checked$LLRparent2)) # no na's

library(dplyr)
library(tidyr)
library(ggplot2)

trios_test_long <- trios_test_checked %>%
  pivot_longer(cols = c(LLRparent1, LLRparent2),
               names_to = "ParentNum",
               values_to = "LLR"
               )
head(trios_test_long)

#################################################################################
# Stats and Plots for GMR Test
# Assign: 79/89 (88.76%) Accuracy: 89/95 (93.68%) Composite: 83.158
ggplot(trios_test_long, aes(x = valid_cross, y = LLR, fill = valid_cross)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.3) +
  labs(
    x = "Were inferred parents truly crossed?",
    y = "LLR(Parent)",
    title = "Distribution of Individual Parent LLR Scores (GMR Test Group)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )

ggplot(trios_test_checked, aes(x = valid_cross, y = LLRpair, fill = valid_cross)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.3) +
  labs(
    x = "Were inferred parents truly crossed?",
    y = "LLR(Pair)",
    title = "Distribution of Parent Pair LLR Scores (GMR Test Group)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )
#################################################################################

# to do next: 
# validate the crosses from seq_test and plot LLR's. Compare to gmr_test
# see module and complex test above. they're not super different. goal here is to
# get a bunch of assignments for the test group, USING the settings we used for the
# gmr test, and see if sequoia and GMR are different. 

# THEN, we want to take the module and complex settings and change them to what
# we'll use for the whole dataset, see how THAT differs...

# and FINALLY, then take those settings, apply them to the whole dataset, Module = "ped"
# and complex = "full", and we want to see how those stack up to what we set for gmr_test,
# since that's what we know should have made that test situation as realistic as possible.
# Module = "par" since we know those first order relationships should dominate, and 
# complex = "simp" since we know that to be a monogamous mating structure. 
# i.e., here's what scores should look like in an ideal scenario. when we expand
# to the full set, here's what they look like, but keep in mind we did have to change
# some of these parameters to accomodate for the different biological situation.


#################################################################################

##### ----- ---- #####

##### ----- Notes ---- #####

# notes:
# work on generation times and making the ageprior more informative
# saying you get FS at agedif = 2. ask Jisca? just try some more stuff out first.
  # not a huge impact when discrete gens are defined. see what else needs to be done.

# could try subsetting to just F0 and F1, seeing if we can extract more relationships,
# similar to what we were doing with the test set? do pairwise generations?
  # could ask jisca her thoughts.


# duplicate check? see the seq objects
  # (nothing to sweat. just some inds that are similar) would expect them to be siblings or something. this is a low div pop.
  
  
# per-sample md? DONE
  # 6 indivs scored for < 20 % of snps
  # 2 indivs genotyped for < 5 % of snps
  # should be fixed now.

