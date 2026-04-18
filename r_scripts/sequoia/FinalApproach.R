##### ---- FinalApproach.R created by SPJ 022426, edited 041726 #####
## PURPOSE: over the past few weeks, i have come up with this two-pronged approach 
# that relies on 1. CalcPairLL() -> LLtoProb() and 2. sequoia() -> GetMaybeRel()
# -> GetRelM() to 1. come up with probabilities that each pair of inds is related
# according to each possible relationship type, and 2. infer the most likely pedigree
# and generate a relatedness matrix from that. (see extrapolation_120925.R)

# because it'll be quite difficult to articulate that approach and to interpret
# results in a talk or paper format, we've decide it's best to combine the two approaches,
# or to pick one of the two. it appears as though you can run CalcPairLL() conditional
# on an inferred pedigree, so the plan is to do that on each of the groups and then
# make the same plots as before. number of PO assignments. 

# 04/17/26
# okay, so THAT^ approach did not work. the pedigree that you use as your "conditional"
# one needs to be from real, observational data. so combining the two approahces
# (doing pairwise and conditioning on the pedigree output from sequoia() and GetMaybeRel())
# wasn't feasible. we went with the approach in SaugerParentage/r_scripts/sequoia/
# extrapolation_120925.R, and produced some figures from that, but we're no longer
# satisfied with those, and we're looking for something more informative. in this
# script, i'm starting with a CLEAN environment, and attempting to run things from
# scratch to produce a stacked barplot and a big boxplot with all PO and GP stuff
# on there. we start with the genotype matrix for all individuals.

setwd("/Users/samjohnson/Documents/Sauger_102325/GeneticData/Sequoia/Sequoia_Inp/contam_fastp_svit_mem/firstfilt_hardfilt_thinned/mindep8_maf30/geno_mat")
getwd()

##### ---- packages #####
library(sequoia)
library(tidyverse)
library(ggplot2)
library(ggpattern)

##### ---- load and initial prep for genotype matrix #####

mat_thin100K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012_conv",
                           header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

# remove first column. ascending numbers 1 - 1184
mat_thin100K <- mat_thin100K[, -1]

# convert to matrix
mat_thin100K <- as.matrix(mat_thin100K)

dim(mat_thin100K)

# read in the sample id's
ind100K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012.indv", header = FALSE)

# rename first column to something meaningful
ind100K <-  ind100K %>% 
  rename(sample = V1)

# establish sample identities in the geno_mat
rownames(mat_thin100K) <- ind100K$sample

# read in scaffolds and positions
pos100K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012.pos", header = FALSE)

# create full positions by combining the two columns
pos100K$position <-  paste(pos100K$V1, pos100K$V2, sep = "_")

# set positions as column names of geno_mat
colnames(mat_thin100K) <- pos100K$position

##### ----- ---- #####

##### ---- adding Additional Data, Filtering Individuals (important!!) #### 

# --- LH FOR ALL --- #
# first i need to add LH data, so that I can filter the GenoM to these indivs.

# IMPORTANT:  UNTIL NOW, I HAVE DONE THIS WITH SOME F0 INDIVIDUALS THAT WERE NOT
# CROSSED. ADDITIONALLY, I WANT TO REMOVE THE JUVENILE F1'S. GOING TO DO BOTH OF
# THOSE STEPS HERE.

LH_All <- read.csv(file = "LH_F0_F1Spawn_F1Juv_F2.csv", header = TRUE)
# change the F's to 1 and M's to 2, all others are 3's
LH_All$Sex[LH_All$Sex == "M"] <- 2
LH_All$Sex[LH_All$Sex == "F"] <- 1
LH_All$Sex <- as.numeric(LH_All$Sex)
str(LH_All)
dim(LH_All)

# FILTERING FOR F0 INDS
LH_F0 <- LH_All %>% 
    filter(BirthYear == 0) # 263, all inds. need to filter down to just
                           # the ones that were crossed.

sar_2015_filt_split_pair <- read.csv(file = "sar_2015_filt_split_pair.csv") # 122
# 53F + 59M, arranged into 61 unique crosses ( x2 is 122)
sar_2016_filt_split_pair <- read.csv(file = "sar_2016_filt_split_pair.csv") # 112
# 47F + 55M, arranged into 56 unique crosses ( x2 is 112)

# checks out with info from Documents/Sauger_102325/GeneticData/F0_CROSSES_021326/unique_F0_crosses.R

# okay great so now we filter for only those.
f0s_that_were_crossed <- c(sar_2015_filt_split_pair$Sample_ID, sar_2016_filt_split_pair$Sample_ID) #234
f0s_that_were_crossed <- unique(f0s_that_were_crossed) # 214 HOLY SHIT YES
f0s_that_were_crossed <- data.frame(f0s_that_were_crossed)

f0s_2015 <- f0s_that_were_crossed[1:112,]
f0s_2016 <- f0s_that_were_crossed[113:nrow(f0s_that_were_crossed),]

LH_F0_true <- LH_All %>% 
  filter(ID %in% f0s_that_were_crossed$f0s_that_were_crossed) # 214

# FILTERING FOR F1 INDS
LH_F1 <- LH_All %>% 
  filter(BirthYear == 1) %>%  # 334
  filter(!grepl("SAR_15_", ID)) # 315, checks out w/ SAR_Data_021026

LH_F2 <- LH_All %>% 
  filter(BirthYear == 2) #480 # checks out w/ SAR_Data_021026

LH_All <- bind_rows(LH_F0_true, LH_F1, LH_F2) # 1009 inds
table(LH_All$Sex)
# 114 males, 100 females from both years.
# 53 + 47 females = 100
# 59 + 55 males = 114

# --- --- # 

table(LH_All$BirthYear)
# 214 F0's (15 and 16)
# 315 F1's (spawning adults from 2021)
# NO TEST F1's included
# 480 F2's (19, 20, 21, all over)
# = 1009
# note: will's samples excluded, possible hybrids excluded, WF fish excluded




# --- LH FOR TEST --- #
# read in parent offspring (F0 and Test F1 sample id's)
testsamp <- read.csv(file = "posampleids.csv", header = TRUE) # 210 samples, 114 Parents, 96 Test F1's

# read in lh data for test inds
LH_Test <- read.csv(file = "testindivs_LH.csv", header = TRUE)
LH_Test_Parents <- LH_Test %>% 
    filter(ID %in% f0s_2015) # 112
LH_Test_Offspring <- LH_Test[1:95,]
LH_Test <- bind_rows(LH_Test_Parents, LH_Test_Offspring) # 207 inds. great.

# change the F's to 1 and M's to 2
LH_Test$Sex[LH_Test$Sex == "M"] <- 2
LH_Test$Sex[LH_Test$Sex == "F"] <- 1
LH_Test$Sex <- as.numeric(LH_Test$Sex)
str(LH_Test)
table(LH_Test$Sex) # 53 female F0, 59 male F0, 95 unknown sex test F1
LH_Test$BirthYear[1:112] <- 0
LH_Test$BirthYear[113:nrow(LH_Test)] <- 1
# THIS IS READY NOW FOR GETMAYBEREL



# --- Filtering the mat_thin100K GenoM to include each of the groups of inds --- #

# filter GenoM so it only includes ids from LH_All$ID
dim(mat_thin100K) # 1184, 943
mat_thin100K_all <- mat_thin100K[rownames(mat_thin100K) %in% LH_All$ID, , drop = FALSE]
dim(mat_thin100K_all) # 992 943
dim(LH_All)       # 1009 3

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
dim(mat_thin100K_all) # 992 943
dim(LH_All)       # 992 3
# okay, that's remedied. great.

# filter GenoM so it only includes ids from LH_Test$ID
dim(mat_thin100K) # 1184, 943
mat_thin100K_test <- mat_thin100K[rownames(mat_thin100K) %in% LH_Test$ID, , drop = FALSE]
dim(mat_thin100K_test) # 207 943
dim(LH_Test)       # 207 3

# so now we're dealing with two pairs of objects: mat_thin100K_all with LH_All
# and mat_thin100K_test with LH_Test

dim(mat_thin100K_all)
dim(LH_All)

dim(mat_thin100K_test)
dim(LH_Test)

# nice.
##### ----- ---- #####

##### ---- HWE Filtering and Final GenoMat Maintenence #####
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
# ℹ After exclusion, There are  990  individuals and  943  SNPs.



# WORK ON TEST DATASET (mat_thin100K_test)
# check for all heterozygous sites since those likely will not be informative.
all_ones_100_test <- apply(mat_thin100K_test, 2, function(col) all(col == 1))
# remove all het. sites
mat_thin100K_test <- mat_thin100K_test[, !all_ones_100_test]
dim(mat_thin100K_test) # lost one

check_thin100K_test <- CheckGeno(mat_thin100K_test, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
# ✔ Genotype matrix looks OK! There are  207  individuals and  942  SNPs.
##### ----- ---- #####

##### ---- Missing Data per Individual #####
# first pass says that there are six individuals genotyped for <20% of snps. 
# let's remedy that.

miss_per_ind <- read.table(file = "miss_per_indv_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.imiss", header = TRUE)
dim(miss_per_ind)
inds_to_keep <- miss_per_ind$INDV[miss_per_ind$F_MISS < 0.2]
inds_to_discard <- miss_per_ind$INDV[miss_per_ind$F_MISS > 0.2]

# individuals who are genotyped for 80% or samples or MORE
length(inds_to_keep) # 1154
dim(miss_per_ind) # 1184, means 30 inds are filtered out...

dim(check_thin100K_all) # 990 samples from JUST the groups we're interested in (F0, F1Spawn, F2)
check_thin100K_all <- check_thin100K_all[rownames(check_thin100K_all) %in% inds_to_keep, ]
dim(check_thin100K_all) # 969 samples remaining.

dim(check_thin100K_test) # 207 inds
check_thin100K_test <- check_thin100K_test[rownames(check_thin100K_test) %in% inds_to_keep, ]
dim(check_thin100K_test) # lost two inds... somehow...

dim(check_thin100K_all) # 969 inds
dim(check_thin100K_test) # 205 inds

##### ---- ---- #####

# see /Users/samjohnson/Documents/Sauger_102325/GeneticData/Sequoia/...
# sequoia_params/radseq_errsequoia_params_radseq_FINAL.xlsx for data on
# optimizing the rad seq error parameters. here's what i landed on.
errM <- Err_RADseq(E0 = 0.075, E1 = 0.025, Return = 'matrix')

##### ---- sequoia() -> GetMaybeRel() #####
dim(check_thin100K_test) # 206 942
dim(LH_Test) # 208 3
LH_Test <- LH_Test %>% 
  filter(ID %in% rownames(check_thin100K_test))
dim(LH_Test) # 206 3, needed to account for that 
table(LH_Test$BirthYear) # 111 parents, 95 test offspring

seq_test <- sequoia(GenoM = check_thin100K_test,
                    LifeHistData = LH_Test,
                    Module = 'ped',
                    Err = errM,
                    Complex = "full",
                    Herm = "no",
                    UseAge = "yes",
                    args.AP=list(Discrete = TRUE, 
                                 MinAgeParent = 1, MaxAgeParent = 1), # notice here! won't assign F0-F2 PO
                    CalcLLR = TRUE,
                    StrictGenoCheck = TRUE,
                    DummyPrefix = c("F", "M"),
                    Tfilter = -2,
                    Tassign = 0.5)
# ✔ assigned 98 dams and 98 sires to 206 + 6 individuals (real + dummy) 

gmr_test <- GetMaybeRel(GenoM = check_thin100K_test,
                        SeqList = seq_test,
                        AgePrior = seq_test[["AgePriors"]],
                        Err = errM,
                        Module = "ped",
                        Complex = "full",
                        Herm = "no",
                        quiet = FALSE,
                        Tfilter = -2,
                        Tassign = 0.5,
                        MaxPairs = 20 * nrow(check_thin100K_test))
# ✔ Found 0 likely parent-offspring pairs, and 42, other non-assigned pairs of possible relatives

seq_test_summary <- SummarySeq(SeqList = seq_test)

cpll_test <- CalcPairLL(GenoM = check_thin100K_test,
                        LifeHistData = LH_Test,
                        AgePrior = TRUE,
                        SeqList = seq_test,
                        Module = "ped",
                        Complex = "full",
                        Herm = "no",
                        quiet = FALSE,
                        InclDup = FALSE,
                        Err = errM,
                        Tassign = 0.5,
                        Tfilter = -2,
                        Plot = TRUE)

IDs <- rownames(check_thin100K_test)
Pairs_test <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_test <- Pairs_test %>% 
  left_join(LH_Test %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_Test %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_test$AgeDif <- Pairs_test$BY2 - Pairs_test$BY1

Pairs_test$focal <- "U"

dim(Pairs_test)

PairLL_test <- CalcPairLL(Pairs = Pairs_test,
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
prob_pairs_test <- plyr::aaply(as.matrix(PairLL_test[,10:16]), .margin = 1, LLtoProb)
prob_pairs_test <- cbind(PairLL_test[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_test)



##### ---- ---- #####

