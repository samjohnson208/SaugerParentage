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

# FILTERING FOR F0 INDS, keep only those that were crossed
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

# FILTERING FOR F1 INDS, remove juvenile F1's sampled in fall of 2015
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
LH_All <- LH_All %>% 
  filter(ID %in% rownames(check_thin100K_all))
dim(LH_All) # 969
LH_Test <- LH_Test %>% 
  filter(ID %in% rownames(check_thin100K_test))
dim(LH_Test) # 205

##### ---- ---- #####

# see /Users/samjohnson/Documents/Sauger_102325/GeneticData/Sequoia/...
# sequoia_params/radseq_err/sequoia_params_radseq_FINAL.xlsx for data on
# optimizing the rad seq error parameters. here's what i landed on.
errM <- Err_RADseq(E0 = 0.075, E1 = 0.025, Return = 'matrix')

# STARTING DATAFRAMES
dim(check_thin100K_all) # 969 inds, 943 snps
dim(check_thin100K_test) # 205 inds, 942 snps
dim(LH_All) # 969 inds
dim(LH_Test) # 205 inds

##### ---- make filtered LH dataframes ---- #####
LH_f0 <- LH_All %>% 
  filter(BirthYear == 0)
dim(LH_f0) # 206

LH_f1 <- LH_All %>% 
  filter(BirthYear == 1)
dim(LH_f1) # 309

LH_f2 <- LH_All %>% 
  filter(BirthYear == 2)
dim(LH_f2) # 454

LH_f0f1 <- rbind(LH_f0, LH_f1)
dim(LH_f0f1) # 515

LH_f1f2 <- rbind(LH_f1, LH_f2)
dim(LH_f1f2) # 763

LH_f0f2 <- rbind(LH_f0, LH_f2)
dim(LH_f0f2) # 660

##### ---- make filtered genotype matrices ---- #####
f0_inds <- LH_f0$ID # 206
f1_inds <- LH_f1$ID # 309
f2_inds <- LH_f2$ID # 454
f0f1_inds <- LH_f0f1$ID # 515
f1f2_inds <- LH_f1f2$ID # 763
f0f2_inds <- LH_f0f2$ID # 660

# can't use filter because it's a matrix and not a dataframe...
check_thin100K_f0 <- check_thin100K_all[rownames(check_thin100K_all) %in% f0_inds,] # 206 inds
check_thin100K_f1 <- check_thin100K_all[rownames(check_thin100K_all) %in% f1_inds,] # 309 inds
check_thin100K_f2 <- check_thin100K_all[rownames(check_thin100K_all) %in% f2_inds,] # 454 inds

check_thin100K_f0f1 <- check_thin100K_all[rownames(check_thin100K_all) %in% f0f1_inds,] # 515 inds
check_thin100K_f1f2 <- check_thin100K_all[rownames(check_thin100K_all) %in% f1f2_inds,] # 763 inds
check_thin100K_f0f2 <- check_thin100K_all[rownames(check_thin100K_all) %in% f0f2_inds,] # 660 inds
##### ---- ---- #####

# THIS CODE PULLED FROM SaugerParentage/r_scripts/sequoia/extrapolation_120925.R
# going to rerun everything that produces the dataframes with probabilities. (i.e.,
# the output from the pairwise method). notice how i've commented out the use of 
# GetRelM(), PlotRelPairs(), and SummarySeq(). then we need to find a way to convert 
# those into the plots that we've drawn up. (see photos from 0417/041826)

##### ---- run sequoia(), GetMaybeRel(), GetRelM(), PlotRelPairs() ---- #####
##### ---- f0_f0 ---- #####
seq_f0 <- sequoia(GenoM = check_thin100K_f0,
                  LifeHistData = LH_f0,
                  Module = "ped",
                  Err = errM,
                  Complex = "full",
                  Herm = "no",
                  UseAge = "yes",
                  args.AP=list(Discrete = TRUE, 
                               MinAgeParent = 1, MaxAgeParent = 1),
                  CalcLLR = TRUE,
                  StrictGenoCheck = TRUE,
                  DummyPrefix = c("F", "M"),
                  Tfilter = -2,
                  Tassign = 0.5)
# assigned 24 dams and 24 sires to 206 + 24 individuals (real + dummy)

gmr_f0 <- GetMaybeRel(GenoM = check_thin100K_f0,
                      SeqList = seq_f0,
                      AgePrior = seq_f0[["AgePriors"]],
                      Err = errM,
                      Module = "ped",
                      Complex = "full",
                      LifeHistData = LH_f0,
                      Herm = "no",
                      quiet = FALSE,
                      Tfilter = -2,
                      Tassign = 0.5,
                      MaxPairs = 7 * nrow(check_thin100K_f0))
# Found 0 likely parent-offspring pairs, and 77, other non-assigned pairs of possible relatives

# relm_f0 <- GetRelM(Pedigree = seq_f0[["Pedigree"]],
#                    Pairs = gmr_f0$MaybeRel,
#                    GenBack = 1, 
#                    patmat = FALSE,
#                    directed = TRUE,
#                    Return = 'Matrix')
# table(unique(relm_f0))
# 
# relmf0_plot <- PlotRelPairs(RelM = relm_f0, 
#                             drop.U = TRUE, 
#                             pch.symbols = TRUE,
#                             cex.axis = 0.3,
#                             mar = c(5, 5, 1, 8))
# 
# seq_f0_summary <- SummarySeq(SeqList = seq_f0)
##### ---- f0_f1 ---- #####
seq_f0f1 <- sequoia(GenoM = check_thin100K_f0f1,
                    LifeHistData = LH_f0f1,
                    Module = "ped",
                    Err = errM,
                    Complex = "full",
                    Herm = "no",
                    UseAge = "yes",
                    args.AP=list(Discrete = TRUE, 
                                 MinAgeParent = 1, MaxAgeParent = 1),
                    CalcLLR = TRUE,
                    StrictGenoCheck = TRUE,
                    DummyPrefix = c("F", "M"),
                    Tfilter = -2,
                    Tassign = 0.5)
# assigned 65 dams and 69 sires to 515 + 38 individuals (real + dummy)

gmr_f0f1 <- GetMaybeRel(GenoM = check_thin100K_f0f1,
                        SeqList = seq_f0f1,
                        AgePrior = seq_f0f1[["AgePriors"]],
                        Err = errM,
                        Module = "ped",
                        Complex = "full",
                        LifeHistData = LH_f0f1,
                        Herm = "no",
                        quiet = FALSE,
                        Tfilter = -2,
                        Tassign = 0.5,
                        MaxPairs = 7 * nrow(check_thin100K_f0f1))
# Found 1 likely parent-offspring pairs, and 327, other non-assigned pairs of possible relatives

# relm_f0f1 <- GetRelM(Pedigree = seq_f0f1[["Pedigree"]],
#                      Pairs = gmr_f0f1$MaybeRel,
#                      GenBack = 1, 
#                      patmat = FALSE,
#                      directed = TRUE,
#                      Return = 'Matrix')
# table(unique(relm_f0f1))
# 
# relmf0f1_plot <- PlotRelPairs(RelM = relm_f0f1, 
#                               drop.U = TRUE, 
#                               pch.symbols = TRUE,
#                               cex.axis = 0.3,
#                               mar = c(5, 5, 1, 8))
# # no way are these plots useful if there are this many inds on each axis.
# 
# seq_f0f1_summary <- SummarySeq(SeqList = seq_f0f1)
 





##### ---- f1_f1 ---- #####

seq_f1 <- sequoia(GenoM = check_thin100K_f1,
                  LifeHistData = LH_f1,
                  Module = "ped",
                  Err = errM,
                  Complex = "full",
                  Herm = "no",
                  UseAge = "yes",
                  args.AP=list(Discrete = TRUE, 
                               MinAgeParent = 1, MaxAgeParent = 1),
                  CalcLLR = TRUE,
                  StrictGenoCheck = TRUE,
                  DummyPrefix = c("F", "M"),
                  Tfilter = -2,
                  Tassign = 0.5)
# assigned 20 dams and 20 sires to 309 + 18 individuals (real + dummy)

gmr_f1 <- GetMaybeRel(GenoM = check_thin100K_f1,
                      SeqList = seq_f1,
                      AgePrior = seq_f1[["AgePriors"]],
                      Err = errM,
                      Module = "ped",
                      Complex = "full",
                      LifeHistData = LH_f1,
                      Herm = "no",
                      quiet = FALSE,
                      Tfilter = -2,
                      Tassign = 0.5,
                      MaxPairs = 7 * nrow(check_thin100K_f1))
# Found 0 likely parent-offspring pairs, and 142, other non-assigned pairs of possible relatives

# relm_f1 <- GetRelM(Pedigree = seq_f1[["Pedigree"]],
#                    Pairs = gmr_f1$MaybeRel,
#                    GenBack = 1, 
#                    patmat = FALSE,
#                    directed = TRUE,
#                    Return = 'Matrix')
# table(unique(relm_f1))
# 
# relmf1_plot <- PlotRelPairs(RelM = relm_f1, 
#                             drop.U = TRUE,
#                             pch.symbols = TRUE,
#                             cex.axis = 0.3,
#                             mar = c(5, 5, 1, 8))
# 
# seq_f1_summary <- SummarySeq(SeqList = seq_f1)


##### ---- f1_f2 ---- #####
seq_f1f2 <- sequoia(GenoM = check_thin100K_f1f2,
                    LifeHistData = LH_f1f2,
                    Module = "ped",
                    Err = errM,
                    Complex = "full",
                    Herm = "no",
                    UseAge = "yes",
                    args.AP=list(Discrete = TRUE, 
                                 MinAgeParent = 1, MaxAgeParent = 1),
                    CalcLLR = TRUE,
                    StrictGenoCheck = TRUE,
                    DummyPrefix = c("F", "M"),
                    Tfilter = -2,
                    Tassign = 0.5)
# assigned 61 dams and 61 sires to 763 + 56 individuals (real + dummy)

gmr_f1f2 <- GetMaybeRel(GenoM = check_thin100K_f1f2,
                        SeqList = seq_f1f2,
                        AgePrior = seq_f1f2[["AgePriors"]],
                        Err = errM,
                        Module = "ped",
                        Complex = "full",
                        LifeHistData = LH_f1f2,
                        Herm = "no",
                        quiet = FALSE,
                        Tfilter = -2,
                        Tassign = 0.5,
                        MaxPairs = 7 * nrow(check_thin100K_f1f2))
# Found 39 likely parent-offspring pairs, and 598, other non-assigned pairs of possible relatives
# Found 1 parent-parent-offspring trios

# relm_f1f2 <- GetRelM(Pedigree = seq_f1f2[["Pedigree"]],
#                      Pairs = gmr_f1f2$MaybeRel,
#                      GenBack = 1, 
#                      patmat = FALSE,
#                      directed = TRUE,
#                      Return = 'Matrix')
# table(unique(relm_f1f2))
# 
# relmf1f2_plot <- PlotRelPairs(RelM = relm_f1f2, 
#                               drop.U = TRUE, 
#                               pch.symbols = TRUE,
#                               cex.axis = 0.3,
#                               mar = c(5, 5, 1, 8))
# 
# seq_f1f2_summary <- SummarySeq(SeqList = seq_f1f2)

##### ---- f2_f2 ---- #####

seq_f2 <- sequoia(GenoM = check_thin100K_f2,
                  LifeHistData = LH_f2,
                  Module = "ped",
                  Err = errM,
                  Complex = "full",
                  Herm = "no",
                  UseAge = "yes",
                  args.AP=list(Discrete = TRUE, 
                               MinAgeParent = 1, MaxAgeParent = 1),
                  CalcLLR = TRUE,
                  StrictGenoCheck = TRUE,
                  DummyPrefix = c("F", "M"),
                  Tfilter = -2,
                  Tassign = 0.5)
# assigned 39 dams and 39 sires to 454 + 36 individuals (real + dummy)


gmr_f2 <- GetMaybeRel(GenoM = check_thin100K_f2,
                      SeqList = seq_f2,
                      AgePrior = seq_f2[["AgePriors"]],
                      Err = errM,
                      Module = "ped",
                      Complex = "full",
                      LifeHistData = LH_f2,
                      Herm = "no",
                      quiet = FALSE,
                      Tfilter = -2,
                      Tassign = 0.5,
                      MaxPairs = 7 * nrow(check_thin100K_f2))
# Found 0 likely parent-offspring pairs, and 253, other non-assigned pairs of possible relatives


# relm_f2 <- GetRelM(Pedigree = seq_f2[["Pedigree"]],
#                    Pairs = gmr_f2$MaybeRel,
#                    GenBack = 1, 
#                    patmat = FALSE,
#                    directed = TRUE,
#                    Return = 'Matrix')
# table(unique(relm_f2))
# 
# relmf2_plot <- PlotRelPairs(RelM = relm_f2, 
#                             drop.U = TRUE,
#                             pch.symbols = TRUE,
#                             cex.axis = 0.3,
#                             mar = c(5, 5, 1, 8))
# 
# seq_f2_summary <- SummarySeq(SeqList = seq_f2)

##### ---- f0_f2 ---- #####
seq_f0f2 <- sequoia(GenoM = check_thin100K_f0f2,
                    LifeHistData = LH_f0f2,
                    Module = "ped",
                    Err = errM,
                    Complex = "full",
                    Herm = "no",
                    UseAge = "yes",
                    args.AP=list(Discrete = TRUE, 
                                 MinAgeParent = 2, MaxAgeParent = 2),
                    CalcLLR = TRUE,
                    StrictGenoCheck = TRUE,
                    DummyPrefix = c("F", "M"),
                    Tfilter = -2,
                    Tassign = 0.5)
# assigned 75 dams and 73 sires to 660 + 61 individuals (real + dummy) 


gmr_f0f2 <- GetMaybeRel(GenoM = check_thin100K_f0f2,
                        SeqList = seq_f0f2,
                        AgePrior = seq_f0f2[["AgePriors"]],
                        Err = errM,
                        Module = "ped",
                        Complex = "full",
                        LifeHistData = LH_f0f2,
                        Herm = "no",
                        quiet = FALSE,
                        Tfilter = -2,
                        Tassign = 0.5,
                        MaxPairs = 7 * nrow(check_thin100K_f0f2))
# Found 0 likely parent-offspring pairs, and 448, other non-assigned pairs of possible relatives

# relm_f0f2 <- GetRelM(Pedigree = seq_f0f2[["Pedigree"]],
#                      Pairs = gmr_f0f2$MaybeRel,
#                      GenBack = 1, 
#                      patmat = FALSE,
#                      directed = TRUE,
#                      Return = 'Matrix')
# table(unique(relm_f0f2))
# 
# relmf0f2_plot <- PlotRelPairs(RelM = relm_f0f2, 
#                               drop.U = TRUE, 
#                               pch.symbols = TRUE,
#                               cex.axis = 0.3,
#                               mar = c(5, 5, 1, 8))
# 
# seq_f0f2_summary <- SummarySeq(SeqList = seq_f0f2)
# # looks like a few PO duos between f0 and f2. crazy.








##### ---- all_gens ---- #####
dim(check_thin100K_all)
dim(LH_All)
table(LH_All$BirthYear)

seq_all <- sequoia(GenoM = check_thin100K_all,
                   LifeHistData = LH_All,
                   Module = "ped",
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
# assigned 106 dams and 108 sires to 969 + 76 individuals (real + dummy) 


gmr_all <- GetMaybeRel(GenoM = check_thin100K_all,
                       SeqList = seq_all,
                       AgePrior = seq_all[["AgePriors"]],
                       Err = errM,
                       Module = "ped",
                       Complex = "full",
                       LifeHistData = LH_All,
                       Herm = "no",
                       quiet = FALSE,
                       Tfilter = -2,
                       Tassign = 0.5,
                       MaxPairs = 20 * nrow(check_thin100K_all)) # needed to increase
# Found 43 likely parent-offspring pairs, and 977, other non-assigned pairs of possible relatives
# Found 1 parent-parent-offspring trios

# relm_all <- GetRelM(Pedigree = seq_all[["Pedigree"]],
#                     Pairs = gmr_all$MaybeRel,
#                     GenBack = 2, 
#                     patmat = FALSE,
#                     directed = FALSE, # CHANGED HERE: all good, will require less cleaning down the road...
#                     Return = 'Matrix')
# table(unique(relm_all))
# 
# relmall_plot <- PlotRelPairs(RelM = relm_all, 
#                              drop.U = TRUE, 
#                              pch.symbols = TRUE,
#                              cex.axis = 0.3,
#                              mar = c(5, 5, 1, 8)) # weird. can't plot it with two gens back...
# 
# seq_all_summary <- SummarySeq(SeqList = seq_all)

##### ---- test ---- #####
seq_test <- sequoia(GenoM = check_thin100K_test,
                    LifeHistData = LH_Test,
                    Module = "ped",
                    Err = errM,
                    Complex = "full",
                    Herm = "no",
                    UseAge = "yes",
                    args.AP=list(Discrete = TRUE, 
                                 MinAgeParent = 1, MaxAgeParent = 1),
                    CalcLLR = TRUE,
                    StrictGenoCheck = TRUE,
                    DummyPrefix = c("F", "M"),
                    Tfilter = -2,
                    Tassign = 0.5)
# assigned 98 dams and 98 sires to 205 + 7 individuals (real + dummy)

gmr_test <- GetMaybeRel(GenoM = check_thin100K_test,                    
                        SeqList = seq_test,
                        AgePrior = seq_test[["AgePriors"]],
                        Err = errM,
                        Module = "ped",
                        Complex = "full",
                        LifeHistData = LH_Test,
                        Herm = "no",
                        quiet = FALSE,
                        Tfilter = -2,
                        Tassign = 0.5,
                        MaxPairs = 7*nrow(check_thin100K_test))
# Found 0 likely parent-offspring pairs, and 42, other non-assigned pairs of possible relatives

# relm_test <- GetRelM(Pedigree = seq_test[["Pedigree"]],
#                      Pairs = gmr_test$MaybeRel,
#                      GenBack = 1, 
#                      patmat = FALSE,
#                      directed = TRUE,
#                      Return = 'Matrix')
# table(unique(relm_test))
# 
# relmtest_plot <- PlotRelPairs(RelM = relm_test,
#                               pch.symbols = TRUE,
#                               mar = c(5, 5, 1, 8))
# 
# seq_test_summary <- SummarySeq(SeqList = seq_test)




##### ---- ---- #####

setwd("/Users/samjohnson/Desktop/FinalApproachResults")
save.image(file = "backup_postseqgmr.RData")

##### ---- Pairs_f0 ---- #####
# creating the pairs df: create all combinations of id1 and 2, remove rows where 
# they're the same
library(tidyr)
IDs <- rownames(check_thin100K_f0)
length(IDs)
Pairs_f0 <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_f0 <- Pairs_f0 %>% 
  left_join(LH_f0 %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_f0 %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_f0$AgeDif <- Pairs_f0$BY2 - Pairs_f0$BY1

Pairs_f0$focal <- "U"

dim(Pairs_f0)
##### ---- Pairs_f0f1 ---- #####
# creating the pairs df: create all combinations of id1 and 2, remove rows where 
# they're the same
IDs <- rownames(check_thin100K_f0f1)
length(IDs)
Pairs_f0f1 <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_f0f1 <- Pairs_f0f1 %>% 
  left_join(LH_f0f1 %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_f0f1 %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_f0f1$AgeDif <- Pairs_f0f1$BY2 - Pairs_f0f1$BY1

Pairs_f0f1$focal <- "U"

dim(Pairs_f0f1)

##### ---- Pairs_f1 ---- #####
IDs <- rownames(check_thin100K_f1)
length(IDs)
Pairs_f1 <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_f1 <- Pairs_f1 %>% 
  left_join(LH_f1 %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_f1 %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_f1$AgeDif <- Pairs_f1$BY2 - Pairs_f1$BY1

Pairs_f1$focal <- "U"

dim(Pairs_f1)

##### ---- Pairs_f1f2 ---- #####
IDs <- rownames(check_thin100K_f1f2)
length(IDs)
Pairs_f1f2 <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_f1f2 <- Pairs_f1f2 %>% 
  left_join(LH_f1f2 %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_f1f2 %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_f1f2$AgeDif <- Pairs_f1f2$BY2 - Pairs_f1f2$BY1

Pairs_f1f2$focal <- "U"

dim(Pairs_f1f2)

##### ---- Pairs_f2---- #####
IDs <- rownames(check_thin100K_f2)
length(IDs)
Pairs_f2 <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_f2 <- Pairs_f2 %>% 
  left_join(LH_f2 %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_f2 %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_f2$AgeDif <- Pairs_f2$BY2 - Pairs_f2$BY1

Pairs_f2$focal <- "U"

dim(Pairs_f2)

##### ---- Pairs_f0f2 ---- #####
IDs <- rownames(check_thin100K_f0f2)
length(IDs)
Pairs_f0f2 <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_f0f2 <- Pairs_f0f2 %>% 
  left_join(LH_f0f2 %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_f0f2 %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_f0f2$AgeDif <- Pairs_f0f2$BY2 - Pairs_f0f2$BY1

Pairs_f0f2$focal <- "U"

dim(Pairs_f0f2)

##### ---- Pairs_test ---- #####
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
##### ---- ---- #####
save.image(file = "backup_prePairs_all.RData")
##### ---- Pairs_all ---- #####
IDs <- rownames(check_thin100K_all)
Pairs_all <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_all <- Pairs_all %>% 
  left_join(LH_All %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_All %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_all$AgeDif <- Pairs_all$BY2 - Pairs_all$BY1

Pairs_all$focal <- "U"

dim(Pairs_all)
##### ---- ---- #####

save.image(file = "backup_preCalcPairLL.RData")

##### ---- Getting LLRs and probs for all relationships ---- #####
##### ---- PairLL_f0 -> prob_pairs_f0  ---- #####

PairLL_f0 <- CalcPairLL(Pairs = Pairs_f0,
                        GenoM = check_thin100K_f0,
                        LifeHistData = LH_f0,
                        AgePrior = seq_f0[["AgePriors"]],
                        Module = "ped",
                        Complex = "full",
                        Herm = 'no',
                        InclDup = FALSE,
                        Err = errM,
                        Tassign = 0.5,
                        Tfilter = -2,
                        quiet = FALSE,
                        Plot = TRUE)
prob_pairs_f0 <- plyr::aaply(as.matrix(PairLL_f0[,10:16]), .margin = 1, LLtoProb)
prob_pairs_f0 <- cbind(PairLL_f0[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_f0)

##### ---- PairLL_f0f1 -> prob_pairs_f0f1 ---- #####

PairLL_f0f1 <- CalcPairLL(Pairs = Pairs_f0f1,
                          GenoM = check_thin100K_f0f1,
                          LifeHistData = LH_f0f1,
                          AgePrior = seq_f0f1[["AgePriors"]],
                          Module = "ped",
                          Complex = "full",
                          Herm = 'no',
                          InclDup = FALSE,
                          Err = errM,
                          Tassign = 0.5,
                          Tfilter = -2,
                          quiet = FALSE,
                          Plot = TRUE)
prob_pairs_f0f1 <- plyr::aaply(as.matrix(PairLL_f0f1[,10:16]), .margin = 1, LLtoProb)
prob_pairs_f0f1 <- cbind(PairLL_f0f1[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_f0f1)

##### ---- PairLL_f1 -> prob_pairs_f1 ---- #####

PairLL_f1 <- CalcPairLL(Pairs = Pairs_f1,
                        GenoM = check_thin100K_f1,
                        LifeHistData = LH_f1,
                        AgePrior = seq_f1[["AgePriors"]],
                        Module = "ped",
                        Complex = "full",
                        Herm = 'no',
                        InclDup = FALSE,
                        Err = errM,
                        Tassign = 0.5,
                        Tfilter = -2,
                        quiet = FALSE,
                        Plot = TRUE)
prob_pairs_f1 <- plyr::aaply(as.matrix(PairLL_f1[,10:16]), .margin = 1, LLtoProb)
prob_pairs_f1 <- cbind(PairLL_f1[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_f1)


##### ---- PairLL_f1f2 -> prob_pairs_f1f2 ---- #####

PairLL_f1f2 <- CalcPairLL(Pairs = Pairs_f1f2,
                          GenoM = check_thin100K_f1f2,
                          LifeHistData = LH_f1f2,
                          AgePrior = seq_f1f2[["AgePriors"]],
                          Module = "ped",
                          Complex = "full",
                          Herm = 'no',
                          InclDup = FALSE,
                          Err = errM,
                          Tassign = 0.5,
                          Tfilter = -2,
                          quiet = FALSE,
                          Plot = TRUE)
prob_pairs_f1f2 <- plyr::aaply(as.matrix(PairLL_f1f2[,10:16]), .margin = 1, LLtoProb)
prob_pairs_f1f2 <- cbind(PairLL_f1f2[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_f1f2)

##### ---- PairLL_f2 -> prob_pairs_f2 ---- #####

PairLL_f2 <- CalcPairLL(Pairs = Pairs_f2,
                        GenoM = check_thin100K_f2,
                        LifeHistData = LH_f2,
                        AgePrior = seq_f2[["AgePriors"]],
                        Module = "ped",
                        Complex = "full",
                        Herm = 'no',
                        InclDup = FALSE,
                        Err = errM,
                        Tassign = 0.5,
                        Tfilter = -2,
                        quiet = FALSE,
                        Plot = TRUE)
prob_pairs_f2 <- plyr::aaply(as.matrix(PairLL_f2[,10:16]), .margin = 1, LLtoProb)
prob_pairs_f2 <- cbind(PairLL_f2[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_f2)

##### ---- PairLL_f0f2 -> prob_pairs_f0f2 ---- #####

PairLL_f0f2 <- CalcPairLL(Pairs = Pairs_f0f2,
                          GenoM = check_thin100K_f0f2,
                          LifeHistData = LH_f0f2,
                          AgePrior = seq_f0f2[["AgePriors"]],
                          Module = "ped",
                          Complex = "full",
                          Herm = 'no',
                          InclDup = FALSE,
                          Err = errM,
                          Tassign = 0.5,
                          Tfilter = -2,
                          quiet = FALSE,
                          Plot = TRUE)
prob_pairs_f0f2 <- plyr::aaply(as.matrix(PairLL_f0f2[,10:16]), .margin = 1, LLtoProb)
prob_pairs_f0f2 <- cbind(PairLL_f0f2[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_f0f2)

##### ---- PairLL_test -> prob_pairs_test ---- #####

# IMPORTANT NOTE: I did the test group AFTER all the rest of the combinations in
# this script, seq_test, gmr_test, and relm_test were all made below in the "a quick 
# check on the thresholds and test group" section

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


##### ---- PairLL_all -> prob_pairs_all ---- #####

PairLL_all <- CalcPairLL(Pairs = Pairs_all,
                         GenoM = check_thin100K_all,
                         LifeHistData = LH_All,
                         AgePrior = seq_all[["AgePriors"]],
                         Module = "ped",
                         Complex = "full",
                         Herm = 'no',
                         InclDup = FALSE,
                         Err = errM,
                         Tassign = 0.5,
                         Tfilter = -2,
                         quiet = FALSE,
                         Plot = TRUE)
# here, you're applying LLtoProb (row by row) to the 10th t0 16th columns of the 
# PairLL_all, (see help for that function, that's where )
prob_pairs_all <- plyr::aaply(as.matrix(PairLL_all[,10:16]), .margin = 1, LLtoProb)
prob_pairs_all <- cbind(PairLL_all[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_all)

save.image(file = "backup_postPairLL_all.RData")

##### ---- after creating those... ---- #####

##### ---- keep only distinct pairs ---- #####

# first, create a "pair" column, and keep only distinct pairs

prob_pairs_f0_unique <- prob_pairs_f0 %>%
  mutate(Pair = paste(pmin(ID1, ID2), pmax(ID1, ID2), sep = "_")) %>%
  distinct(Pair, .keep_all = TRUE)

prob_pairs_f0f1_unique <- prob_pairs_f0f1 %>%
  mutate(Pair = paste(pmin(ID1, ID2), pmax(ID1, ID2), sep = "_")) %>%
  distinct(Pair, .keep_all = TRUE)

prob_pairs_f1_unique <- prob_pairs_f1 %>%
  mutate(Pair = paste(pmin(ID1, ID2), pmax(ID1, ID2), sep = "_")) %>%
  distinct(Pair, .keep_all = TRUE)

prob_pairs_f1f2_unique <- prob_pairs_f1f2 %>%
  mutate(Pair = paste(pmin(ID1, ID2), pmax(ID1, ID2), sep = "_")) %>%
  distinct(Pair, .keep_all = TRUE)

prob_pairs_f2_unique <- prob_pairs_f2 %>%
  mutate(Pair = paste(pmin(ID1, ID2), pmax(ID1, ID2), sep = "_")) %>%
  distinct(Pair, .keep_all = TRUE)

prob_pairs_f0f2_unique <- prob_pairs_f0f2 %>%
  mutate(Pair = paste(pmin(ID1, ID2), pmax(ID1, ID2), sep = "_")) %>%
  distinct(Pair, .keep_all = TRUE)

prob_pairs_test_unique <- prob_pairs_test %>%
  mutate(Pair = paste(pmin(ID1, ID2), pmax(ID1, ID2), sep = "_")) %>%
  distinct(Pair, .keep_all = TRUE)

prob_pairs_all_unique <- prob_pairs_all %>%
  mutate(Pair = paste(pmin(ID1, ID2), pmax(ID1, ID2), sep = "_")) %>%
  distinct(Pair, .keep_all = TRUE)
# OKAY, THE ALL HAS AN ISSUE WITH THE AgeDif. There are -1's in there, but they
# won't affect the GP stuff so let's move on w it for now. 

##### ---- ---- #####

##### ---- assign generations to the pairs ---- #####
assign_gen <- function(x) {
  case_when(
    x %in% f0_inds ~ "F0",
    x %in% f1_inds ~ "F1",
    x %in% f2_inds ~ "F2",
    TRUE ~ NA_character_ # shouldn't be the case, but put an NA if an ind isn't there
  )
}

prob_pairs_f0_unique_gen <- prob_pairs_f0_unique %>%
  mutate(group_ind1 = assign_gen(ID1),
         group_ind2 = assign_gen(ID2),
         # now we'll use our same pmin and pmax setup to create a group_pair col
         group_pair = paste(pmin(group_ind1, group_ind2),
                            pmax(group_ind1, group_ind2),
                            sep = "_"))

prob_pairs_f0f1_unique_gen <- prob_pairs_f0f1_unique %>%
  mutate(group_ind1 = assign_gen(ID1),
         group_ind2 = assign_gen(ID2),
         # now we'll use our same pmin and pmax setup to create a group_pair col
         group_pair = paste(pmin(group_ind1, group_ind2),
                            pmax(group_ind1, group_ind2),
                            sep = "_"))

prob_pairs_f1_unique_gen <- prob_pairs_f1_unique %>%
  mutate(group_ind1 = assign_gen(ID1),
         group_ind2 = assign_gen(ID2),
         # now we'll use our same pmin and pmax setup to create a group_pair col
         group_pair = paste(pmin(group_ind1, group_ind2),
                            pmax(group_ind1, group_ind2),
                            sep = "_"))

prob_pairs_f1f2_unique_gen <- prob_pairs_f1f2_unique %>%
  mutate(group_ind1 = assign_gen(ID1),
         group_ind2 = assign_gen(ID2),
         # now we'll use our same pmin and pmax setup to create a group_pair col
         group_pair = paste(pmin(group_ind1, group_ind2),
                            pmax(group_ind1, group_ind2),
                            sep = "_"))

prob_pairs_f2_unique_gen <- prob_pairs_f2_unique %>%
  mutate(group_ind1 = assign_gen(ID1),
         group_ind2 = assign_gen(ID2),
         # now we'll use our same pmin and pmax setup to create a group_pair col
         group_pair = paste(pmin(group_ind1, group_ind2),
                            pmax(group_ind1, group_ind2),
                            sep = "_"))

prob_pairs_f0f2_unique_gen <- prob_pairs_f0f2_unique %>%
  mutate(group_ind1 = assign_gen(ID1),
         group_ind2 = assign_gen(ID2),
         # now we'll use our same pmin and pmax setup to create a group_pair col
         group_pair = paste(pmin(group_ind1, group_ind2),
                            pmax(group_ind1, group_ind2),
                            sep = "_"))

test_f1_inds <- LH_Test$ID[111:nrow(LH_Test)]
length(test_f1_inds)
assign_gen_test <- function(x) {
  case_when(
    x %in% f0_inds ~ "F0",
    x %in% test_f1_inds ~ "TestF1",
    TRUE ~ NA_character_ # shouldn't be the case, but put an NA if an ind isn't there
  )
}

prob_pairs_test_unique_gen <- prob_pairs_test_unique %>%
  mutate(group_ind1 = assign_gen_test(ID1),
         group_ind2 = assign_gen_test(ID2),
         # now we'll use our same pmin and pmax setup to create a group_pair col
         group_pair = paste(pmin(group_ind1, group_ind2),
                            pmax(group_ind1, group_ind2),
                            sep = "_"))

prob_pairs_all_unique_gen <- prob_pairs_all_unique %>%
  mutate(group_ind1 = assign_gen(ID1),
         group_ind2 = assign_gen(ID2),
         # now we'll use our same pmin and pmax setup to create a group_pair col
         group_pair = paste(pmin(group_ind1, group_ind2),
                            pmax(group_ind1, group_ind2),
                            sep = "_"))


# now we're workin with...
prob_pairs_f0_unique_gen
prob_pairs_f0f1_unique_gen
prob_pairs_f1_unique_gen
prob_pairs_f1f2_unique_gen
prob_pairs_f2_unique_gen
prob_pairs_f0f2_unique_gen
prob_pairs_test_unique_gen
prob_pairs_all_unique_gen

##### ---- ---- #####

##### ---- take just the PO's and the GP's ---- #####
# prob_pairs_f0_unique_gen
# PO_f1 <- prob_pairs_f1_unique_gen
# PO_f2 <- prob_pairs_f2_unique_gen
# PO_f0f2 <- prob_pairs_f0f2_unique_gen

PO_f0f1 <- prob_pairs_f0f1_unique_gen %>%
  filter(TopRel == "PO")
dim(PO_f0f1)

PO_f1f2 <- prob_pairs_f1f2_unique_gen %>% 
  filter(TopRel == "PO")
dim(PO_f1f2)

PO_test <- prob_pairs_test_unique_gen %>% 
  filter(TopRel == "PO")
dim(PO_test)

GP_all <- prob_pairs_all_unique_gen %>% 
  filter(TopRel == "GP")
dim(GP_all)

##### ---- ---- #####

#### ---- how many assignments and unique individuals in each gen? ---- ####

# here are the four objects that are going to constitute the plots.
PO_test
PO_f0f1
PO_f1f2
GP_all

### PO_test ### 
dim(PO_test) # 172 assignments made
length(unique(PO_test$ID2)) # 94 test f1's assigned
table(table(PO_test$ID2)) # 78 assigned to two, 16 assigned to one

### PO_f0f1 ##
dim(PO_f0f1) # 53 assignments made
length(unique(PO_f0f1$ID2)) # 39 spawning f1's assigned
table(table(PO_f0f1$ID2)) # 1 assigned to three, 12 assigned to two, 26 assigned to one

### PO_f1f2 ##
dim(PO_f1f2) # 39 assignments made
length(unique(PO_f1f2$ID2)) # 37 juvenile f2's assigned  
table(table(PO_f1f2$ID2)) # 2 assigned to two, 35 assigned to one

### GP_all ##
dim(GP_all) # 398 assignments made
length(unique(GP_all$ID2)) # 307 juvenile f2's assigned 
table(table(GP_all$ID2)) # 231 assigned to one, 64 to two, 9 to three, 3 to four

##### ---- ---- #####

getwd()
save.image(file = "backup_prevalid.RData")

# got through recreating all of the information. now we need to continue to go thru
# these objects to validate the inferred sets of parents for each offspring. going
# to take some time for sure. then we can create (AT LEAST) the stacked barplot
# and the boxplots. gonna be okay.

##### ---- preparing to validate against the known F0 crosses ---- #####

# here are the groups need to subset
# test: 1 and 2 inf. parents.
# f0f1: 1 and 2 inf. parents.
# f1f2: 1 and 2 inf. parents.
# gp: 1, 2, 3, and 4 inf. gp's 

PO_test_counts <- PO_test %>%
  group_by(ID2) %>%
  mutate(n_parents = n()) %>%
  ungroup()
table(PO_test_counts$n_parents) # 1 and 2
PO_test_1 <- PO_test_counts %>% filter(n_parents == 1)
PO_test_2 <- PO_test_counts %>% filter(n_parents == 2)

PO_f0f1_counts <- PO_f0f1 %>%
  group_by(ID2) %>%
  mutate(n_parents = n()) %>%
  ungroup()
table(PO_f0f1_counts$n_parents) # 1, 2, and 3
PO_f0f1_1 <- PO_f0f1_counts %>% filter(n_parents == 1)
PO_f0f1_2 <- PO_f0f1_counts %>% filter(n_parents == 2)
PO_f0f1_3 <- PO_f0f1_counts %>% filter(n_parents == 3)

PO_f1f2_counts <- PO_f1f2 %>%
  group_by(ID2) %>%
  mutate(n_parents = n()) %>%
  ungroup()
table(PO_f1f2_counts$n_parents) # 1 and 2
PO_f1f2_1 <- PO_f1f2_counts %>% filter(n_parents == 1)
PO_f1f2_2 <- PO_f1f2_counts %>% filter(n_parents == 2)

GP_all_counts <- GP_all %>%
  group_by(ID2) %>%
  mutate(n_grandparents = n()) %>%
  ungroup()
table(GP_all_counts$n_grandparents) # 1, 2, 3, and 4
GP_all_1 <- GP_all_counts %>% filter(n_grandparents == 1)
GP_all_2 <- GP_all_counts %>% filter(n_grandparents == 2)
GP_all_3 <- GP_all_counts %>% filter(n_grandparents == 3)
GP_all_4 <- GP_all_counts %>% filter(n_grandparents == 4)
##### ---- ---- #####

##### ---- validate the inferred parent/grandparent crosses ---- #####

# the groups we can validate
# Test: inds assigned to 2 parents
# F0F1: inds assigned to 2 parents
# F1F2: none.
# GP: inds assigned to 2 gp's. check 3 and 4 by hand.
# we can do 3 ad 4 w the check_crosses function!

# load f0 crosses
setwd("/Users/samjohnson/Documents/Sauger_102325/GeneticData/F0_CROSSES_021326")
crosses_2015 <- read.csv(file = "sar_2015_filt_split_pair.csv", header = TRUE)
crosses_2016 <- read.csv(file = "sar_2016_filt_split_pair.csv", header = TRUE)
f0_crosses <- bind_rows(crosses_2015, crosses_2016)

# validate groups of inferred parents
PO_test_2_valid <- PO_test_2 %>%
  group_by(ID2) %>%
  mutate(
    inferred_pair = paste(sort(ID1), collapse = "_"),
    valid_cross = inferred_pair %in% f0_crosses$Pair
  ) %>%
  ungroup()

PO_f0f1_2_valid <- PO_f0f1_2 %>% 
  group_by(ID2) %>%
  mutate(
    inferred_pair = paste(sort(ID1), collapse = "_"),
    valid_cross = inferred_pair %in% f0_crosses$Pair
  ) %>%
  ungroup()

PO_f0f1_3_valid <- PO_f0f1_3 %>%
  group_by(ID2) %>%
  mutate(
    n_parents = n(),
    valid_cross = check_crosses(ID1)
  ) %>%
  ungroup()
table(PO_f0f1_3_valid$valid_cross)

# no validation on f1f2 stuff

GP_all_2_valid <- GP_all_2 %>% 
  group_by(ID2) %>%
  mutate(
    inferred_pair = paste(sort(ID1), collapse = "_"),
    valid_cross = inferred_pair %in% f0_crosses$Pair
  ) %>%
  ungroup()

# check for GP 3 and 4
check_crosses <- function(parents) {
  pairs <- combn(sort(parents), 2, FUN = function(x) paste(x, collapse = "_"))
  any(pairs %in% f0_crosses$Pair)
}

GP_all_3_valid <- GP_all_3 %>%
  group_by(ID2) %>%
  mutate(
    n_parents = n(),
    valid_cross = check_crosses(ID1)
  ) %>%
  ungroup()
table(GP_all_3_valid$valid_cross) # 18 F, 9 T, means 6 inds for F, 3 inds for T

GP_all_4_valid <- GP_all_4 %>%
  group_by(ID2) %>%
  mutate(
    n_parents = n(),
    valid_cross = check_crosses(ID1)
  ) %>%
  ungroup()
table(GP_all_4_valid$valid_cross) # 4 F, 8 T, means 1 ind for F, 2 inds for T

##### ---- ---- #####

##### ---- create proportions for plotting ---- #####

# Test 2T/2F/1/0
table(PO_test_2_valid$valid_cross) # 78 inds T, 0 F
length(unique(PO_test_2_valid$ID1)) # 57 parents
length(unique(PO_test_2_valid$ID2)) # 78 offspring 2T
length(unique(PO_test_1$ID2)) # 16 offspring 1
# proportions:
# 78/95 = 2T
# 0/95 = 2F
# 16/95 = 1
# 1/95 = 0

# F0-F1 2T/2F/1/0
table(PO_f0f1_2_valid$valid_cross) # 12 inds T, 0F
length(unique(PO_f0f1$ID1)) # 41 parents
length(unique(PO_f0f1$ID2)) # 39 offspring
length(unique(PO_f0f1_2_valid$ID2)) # 12 offspring 2T
table(PO_f0f1_2_valid$valid_cross)
length(unique(PO_f0f1_1$ID2)) # 26 offspring 1
length(unique(PO_f0f1_3$ID2)) # 1 offspring 3T
# proportions:
# 1/309 = 3T
# 12/309 = 2T
# 0/309 = 2F
# 26/309 = 1
# 271/309 = 0

# F1-F2 2T/2F/1/0
length(unique(PO_f1f2$ID1)) # 38 parents
length(unique(PO_f1f2$ID2)) # 37 offspring
length(unique(PO_f1f2_2$ID2)) # 2 offspring 2U
length(unique(PO_f1f2_1$ID2)) # 35 offspring 1
# proportions:
# 2/454 = 2U
# 35/454 = 1
# 417/454 = 0

# GP 0/1/2T/2F/3T/3F/
length(unique(GP_all$ID1)) # 75 grandparents
length(unique(GP_all$ID2)) # 307 grandchildren
length(unique(GP_all_1$ID2)) # 231 grandchildren 1
table(GP_all_2_valid$valid_cross) # 58 grandchildren 2F, 6 2T
table(GP_all_3_valid$valid_cross) # 6 grandchildren 3F, 3 3T
# 4F AND 4T NEED TO BE EVALUATED DIFFERENTLY. THIS CODE DOES NOT DO IT CORRECTLY
# need to see if BOTH PAIRS are valid for an individual to be 4T
View(GP_all_4_valid)
# SAR_19_5972
# SAR_20_6060
# SAR_21_6353
# none of these have both inferred GP pairs as valid_cross = TRUE
# proportions:
# 3/454 = 4F
# 6/454 = 3F
# 3/454 = 3T
# 58/454 = 2F
# 6/454 = 2T
# 231/454 = 1
# 147/454 = 0

##### ---- ---- #####

##### ---- create stacked barplot ---- #####
setwd("/Users/samjohnson/Desktop/")
summary_df <- read.csv(file = "stacked_barplot.csv", header = TRUE)
str(summary_df)
summary_df$group <- factor(
  summary_df$group,
  levels = c("Test (PO)", "F0-F1 (PO)", "F1-F2 (PO)", "F0-F2 (GP)"))
summary_df$category <- factor(summary_df$category, 
                              levels = c("0", "1", "2T", "2F", 
                                         "2U", "3T", "3F", "4F"))

color_map <- data.frame(category = c("0", "1", "2T", "2U", "2F", "3T", "3F", "4F"),
                        color = c("#453581", "#31688e", "#4cc26c", "#21a585","darkgreen", 
                        "#f7e225", "#fca537", "#da5a6a"))
color_map <- setNames(color_map$color, color_map$category)
                        
summary_df$category <- factor(
  summary_df$category,
  levels = c("0", "1", "2T", "2U", "2F", "3T", "3F", "4F"))

ggplot(summary_df, aes(x = dummy, y = prop, fill = category)) +
  geom_col(width = 0.8, color = "black") +
  facet_wrap(~group, nrow = 1) +
  scale_fill_manual(
    values = color_map,
    name = "n Parents/Grandparents"
  ) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Proportion of Individuals") +
  xlab(NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank()
  )
##### ---- ---- #####

##### ---- create objects for boxplot ---- #####
# Test (1 and 2T)
PO_test_1_toplot <- PO_test_1 %>% 
    mutate(group = "Test (PO)") %>% 
  mutate(category = "1")
PO_test_2_validT_toplot <- PO_test_2_valid %>% 
    filter(valid_cross == TRUE) %>%
    mutate(group = "Test (PO)") %>% 
    mutate(category = "2T")

# F0-F1 (1, 2T, 3T)
PO_f0f1_1_toplot <- PO_f0f1_1 %>% 
  mutate(group = "F0-F1 (PO)") %>% 
  mutate(category = "1")
PO_f0f1_2_validT_toplot <- PO_f0f1_2_valid %>% 
  filter(valid_cross == TRUE) %>%
  mutate(group = "F0-F1 (PO)") %>% 
  mutate(category = "2T")
PO_f0f1_3_validT_toplot <- PO_f0f1_3_valid %>% 
  filter(valid_cross == TRUE) %>%
  mutate(group = "F0-F1 (PO)") %>% 
  mutate(category = "3T")

# F1-F2 (1, 2U)
PO_f1f2_1_toplot <- PO_f1f2_1 %>% 
  mutate(group = "F1-F2 (PO)") %>% 
  mutate(category = "1")
PO_f1f2_2_toplot <- PO_f1f2_2 %>% 
  mutate(group = "F1-F2 (PO)") %>% 
  mutate(category = "2U")

# F0-F2 (GP)
GP_all_1_toplot <- GP_all_1 %>% 
  mutate(group = "F0-F2 (GP)") %>% 
  mutate(category = "1")
GP_all_2_validT_toplot <- GP_all_2_valid %>% 
  filter(valid_cross == TRUE) %>% 
  mutate(group = "F0-F2 (GP)") %>% 
  mutate(category = "2T")
GP_all_2_validF_toplot <- GP_all_2_valid %>% 
  filter(valid_cross == FALSE) %>% 
  mutate(group = "F0-F2 (GP)") %>% 
  mutate(category = "2F")
GP_all_3_validT_toplot <- GP_all_3_valid %>% 
  filter(valid_cross == TRUE) %>% 
  mutate(group = "F0-F2 (GP)") %>% 
  mutate(category = "3T")
GP_all_3_validF_toplot <- GP_all_3_valid %>% 
  filter(valid_cross == FALSE) %>% 
  mutate(group = "F0-F2 (GP)") %>% 
  mutate(category = "3F")
GP_all_4_validF_toplot <- GP_all_4_valid %>% 
  filter(valid_cross == FALSE) %>% 
  mutate(group = "F0-F2 (GP)") %>% 
  mutate(category = "4F")

boxplot_inp <- bind_rows(PO_test_1_toplot, PO_test_2_validT_toplot, 
                         PO_f0f1_1_toplot, PO_f0f1_2_validT_toplot, PO_f0f1_3_validT_toplot,
                         PO_f1f2_1_toplot, PO_f1f2_2_toplot,
                         GP_all_1_toplot, GP_all_2_validT_toplot, GP_all_2_validF_toplot, 
                         GP_all_3_validT_toplot, GP_all_3_validF_toplot, GP_all_4_validF_toplot)
boxplot_inp <- boxplot_inp %>%
  mutate(
    prob_toprel = case_when(
      TopRel == "PO" ~ PO,
      TopRel == "GP" ~ GP,
      TRUE ~ NA_real_
    )
  )
##### ---- ---- #####

##### ---- create boxplot ---- #####
boxplot_inp$group <- factor(boxplot_inp$group,
                            levels = c("Test (PO)", "F0-F1 (PO)", 
                                       "F1-F2 (PO)", "F0-F2 (GP)"))

boxplot_inp$category <- factor(boxplot_inp$category,
                               levels = c("0", "1", "2T", "2F", 
                                          "2U", "3T", "3F", "4F"))

ggplot(boxplot_inp, aes(x = category, y = prob_toprel, fill = category)) +
  geom_boxplot(width = 0.7, color = "black", outlier.size = 0.8) +
  facet_wrap(~group, nrow = 1, scales = "free_x") +
  scale_fill_manual(
    values = color_map,
    name = "n Parents/Grandparents"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  ylab("Assignment Probability") +
  xlab(NULL) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
##### ---- ---- #####











