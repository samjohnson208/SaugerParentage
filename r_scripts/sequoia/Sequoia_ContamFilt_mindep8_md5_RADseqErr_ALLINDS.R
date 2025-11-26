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

# NOTE: AS OF 111125, ALL OF THIS WORK SO FAR HAS BEEN DONE ON VERSION 3.0.3
# STARTING 111125, I AM INVESTIGATING THE NEWLY RELEASED 3.1.2 WHICH SHOULD HELP
# ACCORDING TO JISCA

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
getwd()

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

seq <- sequoia(GenoM = check_thin100K_all, 
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
seq <- sequoia(GenoM = check_thin100K_90, 
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
samp_loci_90 <- sample(ncol(check_thin100K_all), 849)
check_thin100K_90 <- check_thin100K_all[, samp_loci_90]

samp_loci_80 <- sample(ncol(check_thin100K_all), 754)
check_thin100K_80 <- check_thin100K_all[, samp_loci_80]

samp_loci_70 <- sample(ncol(check_thin100K_all), 660)
check_thin100K_70 <- check_thin100K_all[, samp_loci_70]

samp_loci_60 <- sample(ncol(check_thin100K_all), 566)
check_thin100K_60 <- check_thin100K_all[, samp_loci_60]

samp_loci_50 <- sample(ncol(check_thin100K_all), 472)
check_thin100K_50 <- check_thin100K_all[, samp_loci_50]

samp_loci_40 <- sample(ncol(check_thin100K_all), 377)
check_thin100K_40 <- check_thin100K_all[, samp_loci_40]

samp_loci_30 <- sample(ncol(check_thin100K_all), 283)
check_thin100K_30 <- check_thin100K_all[, samp_loci_30]

samp_loci_20 <- sample(ncol(check_thin100K_all), 189)
check_thin100K_20 <- check_thin100K_all[, samp_loci_20]

samp_loci_10 <- sample(ncol(check_thin100K_all), 94)
check_thin100K_10 <- check_thin100K_all[, samp_loci_10]

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


gmr <- GetMaybeRel(GenoM = check_thin100K_all,
                   Err = errM,
                   Module = "ped",
                   Complex = "full",
                   LifeHistData = LH_All,
                   quiet = TRUE,
                   Tfilter = -2,
                   Tassign = 1.0,
                   MaxPairs = 7 * nrow(check_thin100K_all))

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
LH_All <- LH_All %>% 
  filter(ID %in% rownames(check_thin100K_all))
dim(LH_Test) # 208 inds (same here...)
LH_Test <- LH_Test %>% 
    filter(ID %in% rownames(check_thin100K_test))
# fixed 111525

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
# Module = "ped"
# ✔ assigned 104 dams and 105 sires to 206 + 13 individuals (real + dummy) 

# not incredibly different, which is nice. goal here is to make comparisons.
##### ---- ---- #####


################################################################################
# after this, how different is the sequoia test object than the gmr object? (for par and simp)

# run them both
# investigate seq_test object (write code to investigate these seq objects)
# run the valid_cross pipeline on that gmr object, generate vectors of the LLRs
  # for valid cross pairs and analyze those
################################################################################





##### ---- here's the existing code to investigate/validate crosses from gmr objects ---- #####
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
# assign was  89/95 (93.68%), accuracy was 79/89 (88.76%)
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

# create a dataframe that's pivoted so that we have each parent LLR in a row for 
# each individual (2 rows per unique offspring)
trios_test_long <- trios_test_checked %>%
  pivot_longer(cols = c(LLRparent1, LLRparent2),
               names_to = "ParentNum",
               values_to = "LLR"
               )
head(trios_test_long)
##### ---- ---- #####





######### Can we differentiate true/false assignments in gmr using the LLR's? #############
# Stats and Plots for GMR Test
# Assign: 89/95 (93.68%) Accuracy: 79/89 (88.76%) Composite: 83.158
ggplot(trios_test_long, aes(x = valid_cross, y = LLR, fill = valid_cross)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = 16, outlier.colour = "black", 
               outlier.fill = "black", alpha = 0.3) +
  labs(
    x = "Were inferred parents truly crossed?",
    y = "LLR(Parent)",
    title = "Distribution of Individual Parent LLR Scores (GMR Test Group, Module = Par, Complex = Simp)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )

ggplot(trios_test_checked, aes(x = valid_cross, y = LLRpair, fill = valid_cross)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = 16, outlier.colour = "black", 
                outlier.fill = "black", alpha = 0.3) +
  labs(
    x = "Were inferred parents truly crossed?",
    y = "LLR(Pair)",
    title = "Distribution of Parent Pair LLR Scores (GMR Test Group, Module = Par, Complex = Simp)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )
##### ---- ---- #####





#################################################################################
# To expand to the whole dataset, we need to know if/how the two functions differ
# in either their number of assignments and their calculations of the LLR's...

############ 1. Do n assignments or LLR's differ among functions? ###############
# to do next: 
# 1. validate the crosses from seq_test and plot LLR's. Compare to gmr_test
# see module and complex test above. they're not super different. goal here is to
# get a bunch of assignments for the test group, USING the settings we used for the
# gmr test, and see if sequoia and GMR are different.

gmr_trios_test_checked <- trios_test_checked
# just to make sure that the names are absolutely clear

seq_test <- sequoia(GenoM = check_thin100K_test,
                    LifeHistData = LH_Test,
                    Module = "par",
                    Err = errM,
                    Complex = "simp",
                    Herm = "no",
                    UseAge = "yes",
                    CalcLLR = TRUE,
                    StrictGenoCheck = TRUE,
                    DummyPrefix = c("F", "M"),
                    Tfilter = -2,
                    Tassign = 1.0)
# ✔ assigned 91 dams and 90 sires to 206 individuals 

# so we've got a few cases here with seq_test
# ind is assigned to no parents
# ind is assigned to one dummy parent
# ind is assigned to two dummy parents
# ind is assigned to one real parent
# ind is assigned to one real parent, and one dummy parent
# ind is assigned to two real parents


# first, need to implement valid_cross into this output ^ seq_test

seq_ped_test <- as.data.frame(seq_test[["PedigreePar"]])

# read in lookup table
cross_lookup <- read.csv(file = "true_cross_lookup.csv", header = TRUE)
cross_lookup <- cross_lookup[,-1]

# this is a ridiculously effective set of piped functions here:
seq_ped_test <- seq_ped_test %>%                                            # here (2)
  # create a new column called pair, where the cross is structured the same way
  # as those in the lookup table. again, we're using pmin and pmax. sorting the
  # pair entries this way ensures that the cross A_B will be treated the same as
  # the cross B_A.
  mutate(pair = paste(pmin(dam, sire), pmax(dam, sire), sep = "_")) %>%
  
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
length(unique(seq_ped_test$id[grepl("^SAR_15_67", seq_ped_test$id)]))
# 95. has all 95 test f1's in the table.
table(seq_ped_test$valid_cross)
# 78 true. nice.
table(gmr_test$id)

# need to see if the trios (two real parents) agree with GMR.

colnames(gmr_trios_test_checked)
colnames(seq_ped_test)


# here's a function that takes a df and, for all valid crosses, grabs the distinct trios
standardize_pairs <- function(df, id_col = "id", pair_col = "pair") {
  df %>%
    select(all_of(c(id_col, pair_col, "valid_cross"))) %>%
    filter(valid_cross == TRUE) %>%
    distinct() %>%
    rename(id = all_of(id_col), pair = all_of(pair_col))
}

# then here are the valid trios for each
valid_seq <- standardize_pairs(seq_ped_test)
valid_gmr <- standardize_pairs(gmr_trios_test_checked)
dim(valid_gmr) # 79 3
dim(valid_seq) # 78 3


# overlap gets you a df that's the valid dfs joined by id AND pair, and the nrow
# is the number of trios shared between the two methods
overlap <- inner_join(valid_seq, valid_gmr, by = c("id", "pair"))
n_overlap <- nrow(overlap) # 78. UNREAL. 

# which trios were unique to each method?
unique_seq <- anti_join(valid_seq, valid_gmr, by = c("id", "pair"))
unique_gmr <- anti_join(valid_gmr, valid_seq, by = c("id", "pair"))

# summarize those results. 
tibble(
  Method = c("sequoia()", "GetMaybeRel()"),
  ValidCrosses = c(nrow(valid_seq), nrow(valid_gmr)),
  UniqueCrosses = c(nrow(unique_seq), nrow(unique_gmr)),
  SharedCrosses = n_overlap
)
# GMR picked out one extra valid trio. pretty good. comparable searching of pedigree
# space. does it give comparable LLR's for each pair?

seq_llr <- seq_ped_test %>%
  select(id, pair, valid_cross, LLRpair) %>%
  mutate(method = "sequoia")

gmr_llr <- gmr_trios_test_checked %>%
  select(id, pair, valid_cross, LLRpair) %>%
  mutate(method = "GetMaybeRel")

# combine into one tidy dataframe
llr_combined <- bind_rows(seq_llr, gmr_llr)

head(llr_combined) # it has so many rows because sequoia() doesn't put a -555 placeholder
# in for cases where id is an F0, like GetMaybeRel() does. in the gmr_trios_test_checked
# object, those cases have already been filtered out, but this is not the case for
# the seq_ped_test object. don't sweat this, since those weird cases don't have
# an LLR associated with them anyway, they won't affect the plots.

ggplot(llr_combined, aes(x = valid_cross, y = LLRpair, fill = method)) +
  geom_violin(trim = FALSE, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.15, outlier.shape = 16, outlier.colour = "black",
               outlier.fill = "black", alpha = 0.3, position = position_dodge((width = 0.8))) +
  labs(
    x = "Were inferred parents truly crossed?",
    y = "LLR(Parent)",
    title = "Distribution of Parent Pair LLR Scores sequoia() vs. GetMaybeRel() - Test Group - Module = par, Complex = simp"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )

# CONCLUSIONS: the same crosses (except for 1) are being recovered by gmr and sequoia.
# the LLRs are also the same, and the pattern for differentiation of T and F valid cross
# holds true no matter which function you are using. most of the TRUE valid crosses have
# an LLR pair that is about 18. Most FALSE valid crosses have an LLRpair that's lower, around 11.

##### ---- ---- #####




############ 2. Does this change when we mess with the module or complex arguments? ###############
# 2. THEN, we want to take the module and complex settings and change them to what
# we'll use for the whole dataset, see how THAT differs...

seq_ped_test_CHG <- sequoia(GenoM = check_thin100K_test,
                            LifeHistData = LH_Test,
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
# ✔ assigned 94 dams and 95 sires to 206 + 3 individuals (real + dummy) 

seq_ped_test_CHG <- as.data.frame(seq_ped_test_CHG [["PedigreePar"]])

# read in lookup table
cross_lookup <- read.csv(file = "true_cross_lookup.csv", header = TRUE)
cross_lookup <- cross_lookup[,-1]

# this is a ridiculously effective set of piped functions here:
seq_ped_test_CHG <- seq_ped_test_CHG %>%                              # here (2)
  # create a new column called pair, where the cross is structured the same way
  # as those in the lookup table. again, we're using pmin and pmax. sorting the
  # pair entries this way ensures that the cross A_B will be treated the same as
  # the cross B_A.
  mutate(pair = paste(pmin(dam, sire), pmax(dam, sire), sep = "_")) %>%
  
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
length(unique(seq_ped_test_CHG$id[grepl("^SAR_15_67", seq_ped_test_CHG$id)]))
# 95. has all 95 test f1's in the table.
table(seq_ped_test_CHG$valid_cross)
# 77 true. nice.

# need to see if the trios (two real parents) agree with seq_ped_test

colnames(seq_ped_test_CHG)
colnames(seq_ped_test)


# here's a function that takes a df and, for all valid crosses, grabs the distinct trios
standardize_pairs <- function(df, id_col = "id", pair_col = "pair") {
  df %>%
    select(all_of(c(id_col, pair_col, "valid_cross"))) %>%
    filter(valid_cross == TRUE) %>%
    distinct() %>%
    rename(id = all_of(id_col), pair = all_of(pair_col))
}

# then here are the valid trios for each
valid_seq <- standardize_pairs(seq_ped_test)
valid_seq_CHG <- standardize_pairs(seq_ped_test_CHG)

# overlap gets you a df that's the valid dfs joined by id AND pair, and the nrow
# is the number of trios shared between the two methods
overlap <- inner_join(valid_seq, valid_seq_CHG, by = c("id", "pair"))
n_overlap <- nrow(overlap) # 77. GREAT.

# which trios were unique to each method?
unique_seq <- anti_join(valid_seq, valid_seq_CHG, by = c("id", "pair"))
unique_seq_chg <- anti_join(valid_seq_CHG, valid_seq, by = c("id", "pair"))

# summarize those results. 
tibble(
  Method = c("sequoia(par,simp)", "sequoia(ped,full)"),
  ValidCrosses = c(nrow(valid_seq), nrow(valid_seq_CHG)),
  UniqueCrosses = c(nrow(unique_seq), nrow(unique_seq_chg)),
  SharedCrosses = n_overlap
)
# par,simp picked out one extra valid trio. pretty good. comparable searching of pedigree
# space. does it give comparable LLR's for each pair?

seq_llr <- seq_ped_test %>%
  select(id, pair, valid_cross, LLRpair) %>%
  mutate(method = "parsimp")

seq_chg_llr <- seq_ped_test_CHG %>%
  select(id, pair, valid_cross, LLRpair) %>%
  mutate(method = "pedfull")

# combine into one tidy dataframe
llr_combined_seq <- bind_rows(seq_llr, seq_chg_llr)

head(llr_combined_seq) # it has so many rows because sequoia() doesn't put a -555 placeholder
# in for cases where id is an F0, like GetMaybeRel() does. in the gmr_trios_test_checked
# object, those cases have already been filtered out, but this is not the case for
# the seq_ped_test object. don't sweat this, since those weird cases don't have
# an LLR associated with them anyway, they won't affect the plots.

ggplot(llr_combined_seq, aes(x = valid_cross, y = LLRpair, fill = method)) +
  geom_violin(trim = FALSE, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.15, outlier.shape = 16, outlier.colour = "black",
               outlier.fill = "black", alpha = 0.3, position = position_dodge((width = 0.8))) +
  labs(
    x = "Were inferred parents truly crossed?",
    y = "LLR(Parent)",
    title = "Distribution of Individual Parent LLR Scores sequoia(test group) parsimp vs. pedfull"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )

# CONCLUSIONS: the same crosses (except for 1) are being recovered by parsimp and pedfull
# the LLRs are also the same, and the pattern for differentiation of T and F valid cross
# holds true no matter what the module and breeding system complexity are. most of the 
# TRUE valid crosses have an LLR pair that is about 17 or 18. Most FALSE valid 
# crosses have an LLRpair that's lower, around 11. BUT, the distributions overlap quite
# a bit. the medians are offset, sure, but there is considerable overlap. what happens
# if you get a pairLLR of 15. is it a true or false positive?

##### ---- ---- #####



# SUMMARY SO FAR:
# tests that have been completed so far:
# individual parent LLR scores for test (parsimp)
# parent pair LLR scores for test (parsimp)

# parent pair LLR scores for test sequoia vs gmr (parsimp)
  # i.e., do the number of relationships or the LLR's depend on the function? no.

# individual parent LLR's for test using sequoia() par/simp vs ped/full
  # i.e., do the number of relationships or the LLR's depend on the modules? no.



############# 3. Since we know function and arguments don't matter for the test, what happens when we expand to the whole dataset? #############

# 3. and FINALLY, then take those settings, apply them to the whole dataset, Module = "ped"
# and complex = "full", and we want to see how those stack up to what we set for gmr_test,
# since that's what we know should have made that test situation as realistic as possible.
# Module = "par" since we know those first order relationships should dominate, and 
# complex = "simp" since we know that to be a more simple mating structure than will
# be observed on the actual spawning grounds without human intervention.
# i.e., here's what scores should look like in an ideal scenario. when we expand
# to the full set, here's what they look like, but keep in mind we did have to change
# some of these parameters to accommodate for the different biological situation,
# but that shouldn't be a problem, since even when you do that on the test group,
# it doesn't really change the number of trios that are output, the number of valid
# crosses that are output, or the distributions of the LLR scores.

#  first step is to run sequoia on the full dataset with those ped and full options.

# input object are as follows
dim(check_thin100K_all) # 1030 inds
dim(LH_All) # 1060 inds (have not been filtered for md/ind... does this matter?)

# here's the whole dataset (overlapping generations, pedfull)
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

seq_ped_all <- as.data.frame(seq_all[["PedigreePar"]])

head(seq_ped_all)
# just realized that this is a problem. i'm not sure why, but there are NO assignments
# from any F2 fish to any F1 parents. wonder if this is a problem with the ages,
# so now I'm going to subset to have datasets with all of the combinations and try
# to run sequoia() on just those. 






##### ---- sequoia() on pairwise discrete generations ---- #####
# suppose first i should make the LH datasets for each, and then filter the GenoM
# to include each of those combinations. then use discrete ages for each combination.

inds_all <- read.csv(file = "Inds_F0_F1Spawn_F1Juv_F2.csv", header = TRUE)
table(inds_all$Group)

# F0 and F1
  f0_f1_inds <- inds_all %>% 
      filter(Group %in% c("F0", "F1 - Juvenile", "F1 - Spawning"))
  # 263 f0, 19 f1 juv, 315 f1 spawning, 597 inds
  
  LH_F0_F1 <- LH_All[LH_All$ID %in% f0_f1_inds$ID, , drop = FALSE]
  dim(LH_F0_F1) # 591 inds (going to have gotten rid of some of those missing fin clip ones)
  
  check_thin100K_f0f1 <- check_thin100K_all[rownames(check_thin100K_all) %in% LH_F0_F1$ID, , drop = FALSE]
  dim(check_thin100K_f0f1) # 576 inds. suppose we've lost some of the ones with high md

  seq_f0_f1 <- sequoia(GenoM = check_thin100K_f0f1,
                       LifeHistData = LH_F0_F1,
                       args.AP=list(Discrete = TRUE),
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
# ✔ assigned 69 dams and 66 sires to 576 + 37 individuals (real + dummy) 

  
# F1 and F2
  f1_f2_inds <- inds_all %>% 
      filter(Group %in% c("F1 - Juvenile", "F1 - Spawning", "F2"))
  # 19 f1 juv, 315 f1 spawning, 480 f2, 814 inds
  
  LH_F1_F2 <- LH_All[LH_All$ID %in% f1_f2_inds$ID, , drop = FALSE]
  dim(LH_F1_F2) # 797 inds
  
  check_thin100K_f1f2 <- check_thin100K_all[rownames(check_thin100K_all) %in% LH_F1_F2$ID, , drop = FALSE]
  dim(check_thin100K_f1f2) # 780 inds. suppose we've lost some of the ones with high md
  
  seq_f1_f2 <- sequoia(GenoM = check_thin100K_f1f2,
                       LifeHistData = LH_F1_F2,
                       args.AP=list(Discrete = TRUE, MinAgeParent = 1, MaxAgeParent = 1),
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
# ✔ assigned 53 dams and 53 sires to 780 + 50 individuals (real + dummy) 
# i cannot believe it. there are some sibships, but not a single f2 is assigned
# to a REAL f1 parent.
#
#
  
# F0 and F2
  f0_f2_inds <- inds_all %>% 
    filter(Group %in% c("F0", "F2"))
  # 263 f0, 480 f2, 743 indivs
  
  LH_F0_F2 <- LH_All[LH_All$ID %in% f0_f2_inds$ID, , drop = FALSE]
    dim(LH_F0_F2) # 732 inds
    
  check_thin100K_f0f2 <- check_thin100K_all[rownames(check_thin100K_all) %in% LH_F0_F2$ID, , drop = FALSE]
  dim(check_thin100K_f0f2) #704 inds
  
  seq_f0_f2 <- sequoia(GenoM = check_thin100K_f0f2,
                       LifeHistData = LH_F0_F2,
                       args.AP=list(Discrete = TRUE, MaxAgeParent = 2),
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
# ✔ assigned 61 dams and 61 sires to 704 + 57 individuals (real + dummy)
# as anticipated, this one did not work out. let's try the whole thing w/ discrete gens.
# upped maxageparent to 2
# ✔ assigned 65 dams and 61 sires to 704 + 57 individuals (real + dummy) 

##### ---- ---- #####
  
##### ---- sequoia() on the whole dataset --- #####
# whole dataset (discrete genrations, pedfull)
seq_all <- sequoia(GenoM = check_thin100K_all,
                   LifeHistData = LH_All,
                   args.AP=list(Discrete = TRUE),
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
# ✔ assigned 106 dams and 102 sires to 1030 + 71 individuals (real + dummy) 
# so many questions here. where are the sibships and grandparents and everything?
# are the 333 LLR's placeholders for dummy inds?
# i have MaxAgeParent as 1,1 right now. obviously not biologically meaningful.


seq_ped_all <- as.data.frame(seq_all[["Pedigree"]])
colnames(seq_all)

summary <- SummarySeq(SeqList = seq_all)

# alright, so see the $SibSize data and the plots. There are sibships here. How do
# we check them from the sequoia output if they aren't stored in the seq_all[["Pedigree"]]?
# pretty sure to do that we need to run gmr, conditional on those pedigrees.

# but first, i want to run gmr on those pairwise gens. see if it picks out PO
# duos/trios JUST from the genetic data, since we know that should work

##### ---- ---- ##### 

##### ---- GetMaybeRel() on pairwise gens, w/ and w/o pedigree/ageprior inp ---- #####

# to run: gmr on all inds, w/ and w/o pedigree input, w/ LH, w/ and w/o AgePrior
# gmr on f0_f1, w/ and w/o pedigree input, w/ LH, w/ and w/o AgePrior
# gmr on f1_f2, w/ and w/o pedigree input, w/ LH, w/ and w/o AgePrior
# gmr on f0_f2, w/ and w/o pedigree input, w/ LH, w/ and w/o AgePrior

##### ---- gmr on all inds, w/ and w/o pedigree input, w/ LH, w/ and w/o AgePrior ---- #####

# all inds, with pedigree, with ageprior
gmr_all <- GetMaybeRel(GenoM = check_thin100K_all,
                       SeqList = seq_all,
                       AgePrior = seq_all[["AgePriors"]],
                       Err = errM,
                       Module = "ped",
                       Complex = "full",
                       LifeHistData = LH_All,
                       quiet = FALSE,
                       Tfilter = -2,
                       Tassign = 1.0,
                       MaxPairs = 7 * nrow(check_thin100K_all))
# ℹ Searching for non-assigned relative pairs ... (Module = ped)
# ℹ using Pedigree in SeqList
# ℹ using LifeHist in SeqList
# ℹ using AgePriors in SeqList
# ✔ Genotype matrix looks OK! There are  1030  individuals and  943  SNPs.
# ℹ Conditioning on pedigree with 1101 individuals, 106 dams and 102 sires
# ℹ settings in SeqList$Specs will overrule input parameters
# Transferring input pedigree ...
# Counting opposing homozygous loci between all individuals ...
# Checking for non-assigned relatives ...
# 0   10  20  30  40  50  60  70  80  90  100% 
# |   |   |   |   |   |   |   |   |   |   |
#   ************************WARNING - reached max for maybe-rel, truncated!
#   
#   Checking for Parent-Parent-Offspring trios ...
# ✔ Found 29 likely parent-offspring pairs, and 768, other non-assigned pairs of possible relatives


# all inds, with pedigree, without ageprior
gmr_all <- GetMaybeRel(GenoM = check_thin100K_all,
                       SeqList = seq_all,
                       Err = errM,
                       Module = "ped",
                       Complex = "full",
                       LifeHistData = LH_All,
                       quiet = FALSE,
                       Tfilter = -2,
                       Tassign = 1.0,
                       MaxPairs = 7 * nrow(check_thin100K_all))
# this code just pulls the ageprior from seq_all, so this option won't work unless
# we pull the pedigree from a separately created dataframe. not going to mess w
# that right now.

# all inds, without pedigree, with ageprior
gmr_all <- GetMaybeRel(GenoM = check_thin100K_all,
                       AgePrior = seq_all[["AgePriors"]],
                       Err = errM,
                       Module = "ped",
                       Complex = "full",
                       LifeHistData = LH_All,
                       quiet = FALSE,
                       Tfilter = -2,
                       Tassign = 1.0,
                       MaxPairs = 7 * nrow(check_thin100K_all))
# ℹ Searching for non-assigned relative pairs ... (Module = ped)
# ✔ Genotype matrix looks OK! There are  1030  individuals and  943  SNPs.
# ℹ Not conditioning on any pedigree
# Counting opposing homozygous loci between all individuals ...
# Checking for non-assigned relatives ...
# 0   10  20  30  40  50  60  70  80  90  100% 
# |   |   |   |   |   |   |   |   |   |   |
#   ***********************WARNING - reached max for maybe-rel, truncated!
#   
#   Checking for Parent-Parent-Offspring trios ...
# ✔ Found 89 likely parent-offspring pairs, and 798, other non-assigned pairs of possible relatives
# ✔ Found 14 parent-parent-offspring trios

# did this again. i mean, that's 1060900 pairs.

# all inds, without pedigree, without ageprior
gmr_all <- GetMaybeRel(GenoM = check_thin100K_all,
                       AgePrior = seq_all[["AgePriors"]],
                       Err = errM,
                       Module = "ped",
                       Complex = "full",
                       LifeHistData = LH_All,
                       quiet = FALSE,
                       Tfilter = -2,
                       Tassign = 1.0,
                       MaxPairs = 7 * nrow(check_thin100K_all))
# not going to mess with this. going to have to subset to pairwise generations and
# just go from there.

##### ---- ---- #####


##### ---- gmr on f0_f1, w/ and w/o pedigree input, w/ LH, w/ and w/o AgePrior ---- #####
seq_f0_f1 <- sequoia(GenoM = check_thin100K_f0f1,
                     LifeHistData = LH_F0_F1,
                     args.AP=list(Discrete = TRUE),
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
# ✔ assigned 69 dams and 66 sires to 576 + 37 individuals (real + dummy) 

# gmr_f0_f1 with pedigree, with ageprior
gmr_f0_f1 <- GetMaybeRel(GenoM = check_thin100K_f0f1,
                         SeqList = seq_f0_f1,
                         AgePrior = seq_f0_f1[["AgePriors"]],
                         Err = errM,
                         Module = "ped",
                         Complex = "full",
                         LifeHistData = LH_F0_F1,
                         quiet = FALSE,
                         Tfilter = -2,
                         Tassign = 1.0,
                         MaxPairs = nrow(check_thin100K_f0f1) * nrow(check_thin100K_f0f1))
# (note that i increased the number of pairs that it'll allow)
# ℹ Searching for non-assigned relative pairs ... (Module = ped)
# ℹ using Pedigree in SeqList
# ℹ using LifeHist in SeqList
# ℹ using AgePriors in SeqList
# ✔ Genotype matrix looks OK! There are  576  individuals and  943  SNPs.
# ℹ Conditioning on pedigree with 613 individuals, 69 dams and 66 sires
# ℹ settings in SeqList$Specs will overrule input parameters
# Transferring input pedigree ...
# Counting opposing homozygous loci between all individuals ...
# Checking for non-assigned relatives ...
# 0   10  20  30  40  50  60  70  80  90  100% 
# |   |   |   |   |   |   |   |   |   |   |
#   ****************************************
#   Checking for Parent-Parent-Offspring trios ...
# ✔ Found 4 likely parent-offspring pairs, and 362, other non-assigned pairs of possible relatives

# gmr_f0_f1 with pedigree, without ageprior

# gmr_f0_f1 without pedigree, with ageprior
gmr_f0_f1 <- GetMaybeRel(GenoM = check_thin100K_f0f1,
                         AgePrior = seq_f0_f1[["AgePriors"]],
                         Err = errM,
                         Module = "ped",
                         Complex = "full",
                         LifeHistData = LH_F0_F1,
                         quiet = FALSE,
                         Tfilter = -2,
                         Tassign = 1.0,
                         MaxPairs = nrow(check_thin100K_f0f1) * nrow(check_thin100K_f0f1))
# ℹ Searching for non-assigned relative pairs ... (Module = ped)
# ✔ Genotype matrix looks OK! There are  576  individuals and  943  SNPs.
# ℹ Not conditioning on any pedigree
# Counting opposing homozygous loci between all individuals ...
# Checking for non-assigned relatives ...
# 0   10  20  30  40  50  60  70  80  90  100% 
# |   |   |   |   |   |   |   |   |   |   |
#   ****************************************
#   Checking for Parent-Parent-Offspring trios ...
# ✔ Found 64 likely parent-offspring pairs, and 434, other non-assigned pairs of possible relatives
# ✔ Found 14 parent-parent-offspring trios
# huh. cool.

# gmr_f0_f1 without pedigree, without ageprior
gmr_f0_f1 <- GetMaybeRel(GenoM = check_thin100K_f0f1,
                         Err = errM,
                         Module = "ped",
                         Complex = "full",
                         LifeHistData = LH_F0_F1,
                         quiet = FALSE,
                         Tfilter = -2,
                         Tassign = 1.0,
                         MaxPairs = nrow(check_thin100K_f0f1) * nrow(check_thin100K_f0f1))
# ℹ Searching for non-assigned relative pairs ... (Module = ped)
# ✔ Genotype matrix looks OK! There are  576  individuals and  943  SNPs.
# ℹ Not conditioning on any pedigree
# ℹ Ageprior: Flat 0/1, overlapping generations, MaxAgeParent = 2,2
# Counting opposing homozygous loci between all individuals ...
# Checking for non-assigned relatives ...
# 0   10  20  30  40  50  60  70  80  90  100% 
# |   |   |   |   |   |   |   |   |   |   |
#   ****************************************
#   Checking for Parent-Parent-Offspring trios ...
# ✔ Found 33 likely parent-offspring pairs, and 499, other non-assigned pairs of possible relatives
# ✔ Found 13 parent-parent-offspring trios
# got some in here that are FS with agedif = 1. problem. looks like using the ageprior without
# the pedigree gives you a lot more assignments.
##### ---- ---- #####


##### ---- gmr on f1_f2, w/ and w/o pedigree input, w/ LH, w/ and w/o AgePrior ---- #####
seq_f1_f2 <- sequoia(GenoM = check_thin100K_f1f2,
                     LifeHistData = LH_F1_F2,
                     args.AP=list(Discrete = TRUE),
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
# ✔ assigned 53 dams and 53 sires to 780 + 50 individuals (real + dummy)

# gmr_f1_f2 with pedigree, with ageprior
gmr_f1_f2 <- GetMaybeRel(GenoM = check_thin100K_f1f2,
                         SeqList = seq_f1_f2,
                         AgePrior = seq_f1_f2[["AgePriors"]],
                         Err = errM,
                         Module = "ped",
                         Complex = "full",
                         LifeHistData = LH_F1_F2,
                         quiet = FALSE,
                         Tfilter = -2,
                         Tassign = 1.0,
                         MaxPairs = nrow(check_thin100K_f1f2) * nrow(check_thin100K_f1f2))
# (note that i increased the number of pairs that it'll allow)
# ℹ Searching for non-assigned relative pairs ... (Module = ped)
# ℹ using Pedigree in SeqList
# ℹ using LifeHist in SeqList
# ℹ using AgePriors in SeqList
# ✔ Genotype matrix looks OK! There are  780  individuals and  943  SNPs.
# ℹ Conditioning on pedigree with 830 individuals, 53 dams and 53 sires
# ℹ settings in SeqList$Specs will overrule input parameters
# Transferring input pedigree ...
# Counting opposing homozygous loci between all individuals ...
# Checking for non-assigned relatives ...
# 0   10  20  30  40  50  60  70  80  90  100% 
# |   |   |   |   |   |   |   |   |   |   |
#   ****************************************
#   Checking for Parent-Parent-Offspring trios ...
# ✔ Found 43 likely parent-offspring pairs, and 582, other non-assigned pairs of possible relatives
# ✔ Found 1 parent-parent-offspring trios

# gmr_f1_f2 with pedigree, without ageprior

# gmr_f1_f2 without pedigree, with ageprior
gmr_f1_f2 <- GetMaybeRel(GenoM = check_thin100K_f1f2,
                         AgePrior = seq_f1_f2[["AgePriors"]],
                         Err = errM,
                         Module = "ped",
                         Complex = "full",
                         LifeHistData = LH_F1_F2,
                         quiet = FALSE,
                         Tfilter = -2,
                         Tassign = 1.0,
                         MaxPairs = nrow(check_thin100K_f1f2) * nrow(check_thin100K_f1f2))
# ℹ Searching for non-assigned relative pairs ... (Module = ped)
# ✔ Genotype matrix looks OK! There are  780  individuals and  943  SNPs.
# ℹ Not conditioning on any pedigree
# Counting opposing homozygous loci between all individuals ...
# Checking for non-assigned relatives ...
# 0   10  20  30  40  50  60  70  80  90  100% 
# |   |   |   |   |   |   |   |   |   |   |
#   ****************************************
#   Checking for Parent-Parent-Offspring trios ...
# ✔ Found 44 likely parent-offspring pairs, and 630, other non-assigned pairs of possible relatives
# ✔ Found 1 parent-parent-offspring trios


# gmr_f1_f2 without pedigree, without ageprior
gmr_f1_f2 <- GetMaybeRel(GenoM = check_thin100K_f1f2,
                         Err = errM,
                         Module = "ped",
                         Complex = "full",
                         LifeHistData = LH_F1_F2,
                         quiet = FALSE,
                         Tfilter = -2,
                         Tassign = 1.0,
                         MaxPairs = nrow(check_thin100K_f0f1) * nrow(check_thin100K_f0f1))
# ℹ Searching for non-assigned relative pairs ... (Module = ped)
# ✔ Genotype matrix looks OK! There are  780  individuals and  943  SNPs.
# ℹ Not conditioning on any pedigree
# ℹ Ageprior: Flat 0/1, overlapping generations, MaxAgeParent = 2,2
# Counting opposing homozygous loci between all individuals ...
# Checking for non-assigned relatives ...
# 0   10  20  30  40  50  60  70  80  90  100% 
# |   |   |   |   |   |   |   |   |   |   |
#   ****************************************
#   Checking for Parent-Parent-Offspring trios ...
# ✔ Found 35 likely parent-offspring pairs, and 718, other non-assigned pairs of possible relatives
# ✔ Found 1 parent-parent-offspring trios


# problem with these situations where you use neither is that the ageprior is no longer
# informative and you can have maxageparent set up as for you as something that it shouldn't be.
##### ---- ---- #####


##### ---- gmr on f0_f2, w/ and w/o pedigree input, w/ LH, w/ and w/o AgePrior ---- #####
seq_f0_f2 <- sequoia(GenoM = check_thin100K_f0f2,
                     LifeHistData = LH_F0_F2,
                     args.AP=list(Discrete = TRUE, MaxAgeParent = 2),
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
# ✔ assigned 65 dams and 61 sires to 704 + 57 individuals (real + dummy)  


# gmr_f0_f2 with pedigree, with ageprior
gmr_f0_f2 <- GetMaybeRel(GenoM = check_thin100K_f0f2,
                         SeqList = seq_f0_f2,
                         AgePrior = seq_f0_f2[["AgePriors"]],
                         Err = errM,
                         Module = "ped",
                         Complex = "full",
                         LifeHistData = LH_F0_F2,
                         quiet = FALSE,
                         Tfilter = -2,
                         Tassign = 1.0,
                         MaxPairs = nrow(check_thin100K_f0f2) * nrow(check_thin100K_f0f2))
# (note that i increased the number of pairs that it'll allow)
# ℹ Searching for non-assigned relative pairs ... (Module = ped)
# ℹ using Pedigree in SeqList
# ℹ using LifeHist in SeqList
# ℹ using AgePriors in SeqList
# ✔ Genotype matrix looks OK! There are  704  individuals and  943  SNPs.
# ℹ Conditioning on pedigree with 761 individuals, 61 dams and 61 sires
# ℹ settings in SeqList$Specs will overrule input parameters
# Transferring input pedigree ...
# Counting opposing homozygous loci between all individuals ...
# Checking for non-assigned relatives ...
# 0   10  20  30  40  50  60  70  80  90  100% 
# |   |   |   |   |   |   |   |   |   |   |
#   ****************************************
#   ✔ Found 0 likely parent-offspring pairs, and 531, other non-assigned pairs of possible relatives
# hm....

# gmr_f0_f2 with pedigree, without ageprior

# gmr_f0_f2 without pedigree, with ageprior
gmr_f0_f2 <- GetMaybeRel(GenoM = check_thin100K_f0f2,
                         AgePrior = seq_f0_f2[["AgePriors"]],
                         Err = errM,
                         Module = "ped",
                         Complex = "full",
                         LifeHistData = LH_F0_F2,
                         quiet = FALSE,
                         Tfilter = -2,
                         Tassign = 1.0,
                         MaxPairs = nrow(check_thin100K_f0f2) * nrow(check_thin100K_f0f2))
# ℹ Searching for non-assigned relative pairs ... (Module = ped)
# ✔ Genotype matrix looks OK! There are  704  individuals and  943  SNPs.
# ℹ Not conditioning on any pedigree
# Counting opposing homozygous loci between all individuals ...
# Checking for non-assigned relatives ...
# 0   10  20  30  40  50  60  70  80  90  100% 
# |   |   |   |   |   |   |   |   |   |   |
#   ****************************************
#   ✔ Found 0 likely parent-offspring pairs, and 546, other non-assigned pairs of possible relatives


# gmr_f0_f2 without pedigree, without ageprior
gmr_f0_f2 <- GetMaybeRel(GenoM = check_thin100K_f0f2,
                         Err = errM,
                         Module = "ped",
                         Complex = "full",
                         LifeHistData = LH_F0_F2,
                         quiet = FALSE,
                         Tfilter = -2,
                         Tassign = 1.0,
                         MaxPairs = nrow(check_thin100K_f0f2) * nrow(check_thin100K_f0f2))
# ℹ Searching for non-assigned relative pairs ... (Module = ped)
# ✔ Genotype matrix looks OK! There are  704  individuals and  943  SNPs.
# ℹ Not conditioning on any pedigree
# ℹ Ageprior: Flat 0/1, overlapping generations, MaxAgeParent = 3,3
# Counting opposing homozygous loci between all individuals ...
# Checking for non-assigned relatives ...
# 0   10  20  30  40  50  60  70  80  90  100% 
# |   |   |   |   |   |   |   |   |   |   |
#   ****************************************
#   Checking for Parent-Parent-Offspring trios ...
# ✔ Found 5 likely parent-offspring pairs, and 615, other non-assigned pairs of possible relatives

##### ---- ---- #####


# alright, so we need to understand why we're getting these other relationships with GMR
# is it to supplement? what is the point of using both? it's not like they're returning
# the same relationships in the full datasets between the functions, as they are doing 
# in the test group. is it to supplement the pedigree? how much do the relationships
# overlap across outputs?

# then how can we validate using LLR, given what we found with the test group? given
# the overlap in those distributions of the T/F valid_cross?s

#################################################################################

# work 111125: jisca says we should give the newest version a shot. 3.1.2. install
# sequoia_3.1.2.tar.gz. fixed a few bugs there. also see the new LL2Probs function
# she says interpreting the LLRs is not really straightforward/intuitive.


# work 112425: hey that didn't work at all. binaries are now on cran. 
# got 3.1.3 downloaded and installed. things to do are to check out the 3 scores
# (assignment, accuracy, composite) on the test group and see what error matrix
# performs best. then we look at the assignments for valid cross true and false
# for the tests and we do LLtoProb on those. THEN we go back and expand to the 
# rest of the dataset and see what those probabilities look like!!!

# now we explore the radseq error matrix on this new version.
##### ---- RADseq Error Matrix Exploration ---- #####

# create vectors of e0 and e1 values
e0_vals <- c(0.005, 0.01, 0.025, 0.05, 0.075, 0.1)
e1_vals <- c(0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5)

# create empty dfs to display
assign_rate_mat <- matrix(NA, nrow = length(e1_vals), ncol = length(e0_vals),
                          dimnames = list(paste0("e1_", e1_vals),
                                          paste0("e0_", e0_vals)))
acc_rate_mat <- assign_rate_mat  # same structure

# total number of offspring
n_offspring <- 95

# read in cross lookup table
cross_lookup <- read.csv(file = "true_cross_lookup.csv", header = TRUE)
cross_lookup <- cross_lookup[,-1]

# nested loop for all combinations:
for (i in seq_along(e1_vals)) {
  for (j in seq_along(e0_vals)) {
    
    e1_val <- e1_vals[i]
    e0_val <- e0_vals[j]
    
    # print text to display each combination, track progress
    cat("Running GetMaybeRel for E0 =", e0_val, "and E1 =", e1_val, "\n")
    
    # create the error matrix for that combination
    errM <- Err_RADseq(E0 = e0_val, E1 = e1_val, Return = 'matrix')
    
    # run gmr (change dataset here only! how great is that?!)
    gmr <- GetMaybeRel(GenoM = check_thin100K_test,
                       Err = errM,
                       Module = "par",
                       Complex = "simp",
                       LifeHistData = LH_Test,
                       quiet = TRUE,
                       Tassign = 1.0,
                       MaxPairs = 7 * nrow(check_thin100K_test))
    
    # if there are no trios found, place 0 in that place for assign, NA for accuracy
    if (is.null(gmr[["MaybeTrio"]]) || nrow(gmr[["MaybeTrio"]]) == 0) {
      cat("   No trios found for this combination. Skipping...\n")
      assign_rate_mat[i, j] <- 0
      acc_rate_mat[i, j] <- NA
      next
    }
    
    # how many unique test f1's were assigned to parents?
    assigned_ids <- unique(gmr[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr[["MaybeTrio"]]$id)])
    n_assigned <- length(assigned_ids) # define it as an object
    
    # use that object to calculate assignment rate
    assign_rate <- (n_assigned / n_offspring) * 100 
    
    # make trios an object, create the pair entry, check it against the lookup...
    trios <- gmr[["MaybeTrio"]] %>%
      mutate(pair = paste(pmin(parent1, parent2), pmax(parent1, parent2), sep = "_")) %>%
      left_join(cross_lookup %>%
                  select(pair) %>%
                  distinct() %>%
                  mutate(valid_cross = TRUE),
                by = "pair") %>%
      mutate(valid_cross = ifelse(is.na(valid_cross), FALSE, valid_cross)) %>%
      # ...and filter out any placeholders. 
      filter(!LLRparent1 %in% c(555, -555),
             !LLRparent2 %in% c(555, -555),
             !LLRpair    %in% c(555, -555))
    
    # if we have assigned individuals, grab the number of valid crosses and calculate
    # the accuracy rate by dividing it by the number of unique test f1's that got assigned
    if (n_assigned > 0) {
      valid_ids <- unique(trios$id[trios$valid_cross])
      n_valid <- length(valid_ids)
      acc_rate <- (n_valid / n_assigned) * 100
    } else {
      acc_rate <- NA
    }
    
    # place those values in the correct cell of each summary table.
    assign_rate_mat[i, j] <- assign_rate
    acc_rate_mat[i, j] <- acc_rate
  }
}


# print results, round to three decimal places
cat("\nAssignment rate matrix (%):\n")
print(round(assign_rate_mat, 3))

cat("\nAccuracy rate matrix (%):\n")
print(round(acc_rate_mat, 3))

# create composite score matrix
composite_mat <- (assign_rate_mat * acc_rate_mat) / 100

# use the same row and column names as the other two result tables
dimnames(composite_mat) <- dimnames(assign_rate_mat)

# print results, round to three decimal places
cat("\nComposite score matrix (%):\n")
print(round(composite_mat, 3))

# alright, looks like those same pattens hold. the same error params maximize the
# composite score in the test group, that's E0 = 0.075 and E1 = 0.025
# errM <- Err_RADseq(E0 = 0.075, E1 = 0.025, Return = 'matrix')

##### ---- ---- #####

# work 112425: restored plots and refreshed knowledge on everything i've done already
# looks like this: we've got these overlapping dists of parent LLRs and pair LLRs
# for valid cross T/F for gmr, sequoia, and for both regardless of what arguments
# are used module = par/ped; complex = simp/full. so the LLRs don't tell us anything
# we asked jisca: i get the same numbers of relationships and individuals in the 
# returned relationships from gmr and sequoia. am i using gmr right? what's it for?
# she says they pick up the same, but gmr might not choose one assignment that's
# the most likely for that focal ind. that's what goes to sequoia.
# i also asked how i can interpret the LLRs once i expand to the full dataset.
# how do i tell which crosses are true and false. she pointed me to the new version
# specifically the LLtoProbs function. she also recommended that i retry the test
# group with combinations of E0/E1 for the radseq error matrix. that didn't do anything.
# same combination gave me the same composite score for the test group and this 
# error matrix (check_thin100K_test)

# now that we have that stuff sorted, we want to see what the probabilities are
# from the valid cross T/F ones that come out of gmr for the test group. we have
# data frames established that contain each parent LLR and each pairLLR for the
# trios that come out of gmr. that's this code (above)

    # # create a vector and histogram of the LLRs for each duo and trio
    # head(trios_test_checked)
    # table(is.na(trios_test_checked$LLRparent1)) # no na's
    # table(is.na(trios_test_checked$LLRparent2)) # no na's
    # 
    # library(dplyr)
    # library(tidyr)
    # library(ggplot2)
    # 
    # # create a dataframe that's pivoted so that we have each parent LLR in a row for 
    # # each individual (2 rows per unique offspring)
    # trios_test_long <- trios_test_checked %>%
    #   pivot_longer(cols = c(LLRparent1, LLRparent2),
    #                names_to = "ParentNum",
    #                values_to = "LLR"
    #   )
    # head(trios_test_long)

colnames(trios_test_long)
table(trios_test_long$valid_cross)
# remember, 89 assigned, 79 accurate, so you've got 20 false sets of parent LLR's,
# 10 false pair LLR's, 158 true sets of parent LLR's, and 79 true pair LLR's.

trios_test_long <- trios_test_long %>% 
    mutate(pairprob = LLtoProb(LLRpair)) %>% 
    mutate(parentprob = LLtoProb(LLR))
# only the negative parent LLR's are getting a probability, and the positive parent 
# LLRs and pair LLRs are all NA. why is this?

##### ----- Work 112525 - CalcPairLL() and LLtoProb() ---- #####
# emailed jisca, she says what we're looking for can be achieved by using CalcPairLL()
# and using that output in LLtoProb(). Then you get, for each focal pair of individuals
# the probability that their relationship is each of the five first order relationships.

# first thing to do is to understand what we need for CalcPairLL(), that is:
# the df pairs, (all combinations of test inds, need to figure out how to put sexes in,
# agedif
# genoM
# Pedigree (field pedigree with specified relationships that will be conditioned upon)
# LHData (use normal LH_test)
# age prior (grab from sequoia output)
# SeqList is previous sequoia output, use seqlist$PedigreePar when we want to condition on it.
      # if we do that, it also grabs LHData, AgePriors, and ErrM, overrides input params
      # if not, then we need to specify those elsewhere
# Module = ped
# Complex = full
# Herm = no
# InclDup... eh... idk?
# Err = ErrM
# Tassign = new default for sequoia is 0.5, no longer 1.0
      # (doesn't make a difference for sequoia for test group, i tested it today)
      # still the same number of relationships recovered for this 100K dataset
      # still have the same composite score as normal. don't sweat this.
# Tfilter = -2. Stick w the default. has to be twice as likely to be related than
      # unrelated. think this is fine.
# quiet = FALSE, i want the output messages.
# plot = TRUE

# okay. got the arch. of the command mapped out. question now is, do we give it
# ONLY the combinations of F1-F1, F1-F0 to work with, so we can identify sibships,
# but not make it search all of pedigree space, or do we give it ALL of the combinations,
# and let it do things like identify F0-F0 FS relationships. do we want to just give it
# the combinations to assign PO for F1-F0?

# idk dude. i think we need to give it all of the combinations. let it search pedigree
# space. let it filter out the ones that don't make sense given the age prior. make
# damn sure that we aren't getting FS for F1 to F0. The age prior we set up for the
# test set alread SHOULDN'T let that happen. i'd figure out how to make that pairs df
# and let it fly.

# creating the Pairs df
library(tidyr)
IDs <- rownames(check_thin100K_test)
Pairs <- expand_grid(
  ID1 = IDs,
  ID2 = IDs
) |> 
  dplyr::filter(ID1 != ID2)

Pairs <- Pairs |>
  left_join(LH_Test |> select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) |>
  rename(Sex1 = Sex, BY1 = BirthYear) |>
  
  left_join(LH_Test |> select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) |>
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs$AgeDif <- Pairs$BY2 - Pairs$BY1

Pairs$focal <- "U"

PairLL <- CalcPairLL(Pairs = Pairs,
                     GenoM = check_thin100K_test,
                     LifeHistData = LH_Test,
                     AgePrior = seq_test[["AgePriors"]],
                     Module = "ped",
                     Complex = "full",
                     Herm = 'no',
                     InclDup = FALSE,
                     Err = errM,
                     Tassign = 0.5,
                     Tfilter = -2,
                     quiet = FALSE,
                     Plot = TRUE)
prob_pairs <- plyr::aaply(as.matrix(PairLL[,10:16]), .margin = 1, LLtoProb)
prob_pairs_ids <- cbind(PairLL[, c("ID1", "ID2")], prob_pairs)
head(prob_pairs_ids)

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
  # (nothing to sweat. just some inds that are similar) would expect them to be 
  # siblings or something. this is a low div pop.
  
  
# per-sample md? DONE
  # 6 indivs scored for < 20 % of snps
  # 2 indivs genotyped for < 5 % of snps
  # should be fixed now.









