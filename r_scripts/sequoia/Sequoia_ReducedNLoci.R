## Sequoia_ReducedNLoci.R by SPJ 060225
## PURPOSE: to understand problems assoc. with Sequoia. first step is to attempt
          # running sequoia with reduced numbers of loci. i am concerned that the
          # problems i have been having are due to the use of too many loci given
          # what sequoia specifies in its documentation.
## USAGE: Rscript Sequoia_ReducedNLoci.R

## libraries
################################################################################

# install.packages("sequoia")
# 65
library(sequoia)

# install.packages("dplyr")
# 65
library(dplyr)

################################################################################



## initial steps with reduce n loci
################################################################################
# going to run it locally first. then cluster if run time is too long.
setwd("/Users/samjohnson/Documents/Sauger_042225/GeneticData/genotype_matrices/vcftools--012/")

# read in genotype matrix
mat <- read.table(file = "hard_variants_bial_noindels_q20_mindep8_maxdep75_maf30_miss95.012_conv", header = FALSE, sep = "\t", 
                  na.strings = c("NA", "-1"))

# correct row names as sample id's
mat<- mat[, -1]
gmmat <- as.matrix(mat)
ind <- read.table(file = "hard_variants_bial_noindels_q20_mindep8_maxdep75_maf30_miss95.012.indv", header = FALSE)

str(ind)
ind <- ind %>% 
  rename(sample = V1)
rownames(gmmat) <- ind$sample

# read in parent offspring (F0 and Test F1 sample id's)
testsamp <- read.csv(file = "posampleids.csv", header = TRUE) # 210 samples, 114 Parents, 96 Test F1's

#create genotype matrix parent-offspring by filtering rownames of gmmat to include only those in testsamp$sample
gmmat_po <- gmmat[rownames(gmmat) %in% testsamp$sample, , drop = FALSE]
nrow(gmmat_po) # lost one
unmatched_samples <- testsamp$sample[!testsamp$sample %in% rownames(gmmat)] # SAR_15_6757 wasn't sequenced. (bad dna concentration)
# check genotype matrix for samples/loci to be excluded
gmmat_po_check <- CheckGeno(gmmat_po, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))    
# 95 Test F1's, 114 Parents going into this ^

# run sequoia on the checked genotype matrix
outfull <- sequoia_custom(GenoM = gmmat_po_check, Module = 'ped', StrictGenoCheck = TRUE, CalcLLR = TRUE)
gmr <- GetMaybeRel_Custom(GenoM = gmmat_po_check)

# 24 PO relationships recovered. i won't yet check them against the hiphop retrieved ones
# since 24 is not even close to what we'd need. let's first try subsetting. retain 75% of the loci
# and just ratchet down until we hit 300 or so and see what the numbers do.

# check for all heterozygous sites since those likely will not be informative.
all_ones <- apply(gmmat_po_check, 2, function(col) all(col == 1))
filtered_gmmat_po_check <- gmmat_po_check[, !all_ones]
# removed two sites

# retain 75% of the original loci
samp_loci <- sample(ncol(filtered_gmmat_po_check), 972)
filtered_gmmat_po_check_75 <- filtered_gmmat_po_check[, samp_loci]
outfull75 <- sequoia(GenoM = filtered_gmmat_po_check_75, Module = 'ped', StrictGenoCheck = TRUE, CalcLLR = TRUE)
gmr75 <- GetMaybeRel(GenoM = filtered_gmmat_po_check_75)
# alright, got up to 41 po pairs that time... let's continue...

# retain 50% of the original loci
samp_loci <- sample(ncol(filtered_gmmat_po_check), 648)
filtered_gmmat_po_check_50 <- filtered_gmmat_po_check[, samp_loci]
outfull50 <- sequoia(GenoM = filtered_gmmat_po_check_50, Module = 'ped', StrictGenoCheck = TRUE, CalcLLR = TRUE)
gmr50 <- GetMaybeRel(GenoM = filtered_gmmat_po_check_50)
# got up to 56 po pairs that time...

# retain 40% of the original loci
samp_loci <- sample(ncol(filtered_gmmat_po_check), 520)
filtered_gmmat_po_check_40 <- filtered_gmmat_po_check[, samp_loci]
outfull40 <- sequoia(GenoM = filtered_gmmat_po_check_40, Module = 'ped', StrictGenoCheck = TRUE, CalcLLR = TRUE)
gmr40 <- GetMaybeRel(GenoM = filtered_gmmat_po_check_40)
# also only 56 pairs that time...

# retain 30% of the original loci
samp_loci <- sample(ncol(filtered_gmmat_po_check), 390)
filtered_gmmat_po_check_30 <- filtered_gmmat_po_check[, samp_loci]
outfull30 <- sequoia(GenoM = filtered_gmmat_po_check_30, Module = 'ped', StrictGenoCheck = TRUE, CalcLLR = TRUE)
gmr30 <- GetMaybeRel(GenoM = filtered_gmmat_po_check_30)
# 57 pairs that time...

# next thing to do is to go ahead and check these relationships against what hiphop outputs
# see if the same inferences are being made by both packages.

# would like to check some of these relationships against ones predicted by hiphop
setwd("/Users/samjohnson/Documents/Sauger_042225/GeneticData/hiphop/filtering_exercise/hard_filtered_dataframes/")
test_top_1 <- read.csv(file = "top1_100_75_50_25.csv", header = TRUE)


# maybe it's not returning that many relationships because the min LLR defaults to 0.05 or something?
gmr30 <- GetMaybeRel(GenoM = filtered_gmmat_po_check_30, Tassign = 0)
# nope, still 57 PO's.

# doubt it'll work... but we could try even fewer loci?
samp_loci <- sample(ncol(filtered_gmmat_po_check), 260)
filtered_gmmat_po_check_20 <- filtered_gmmat_po_check[, samp_loci]
outfull20 <- sequoia(GenoM = filtered_gmmat_po_check_20, Module = 'ped', StrictGenoCheck = TRUE, CalcLLR = TRUE)
gmr20 <- GetMaybeRel(GenoM = filtered_gmmat_po_check_20)

# dropped to 260 loci... FINALLY came up w/ 5 trios.
# all five agree with the hiphop output. MAYBE it's the Tfilter argument...
# "more negative values may decrease non-assignment, but will increase computational time"
gmr20 <- GetMaybeRel(GenoM = filtered_gmmat_po_check_20, Tassign = 0, Tfilter = -10)
# 62 likely parent offspring pairs, 20 other pairs, 5 trios.
# idk man... so strange.

# doubt it'll work... but we could try even fewer loci?
samp_loci <- sample(ncol(filtered_gmmat_po_check), 129)
filtered_gmmat_po_check_10 <- filtered_gmmat_po_check[, samp_loci]
outfull10 <- sequoia(GenoM = filtered_gmmat_po_check_10, Module = 'ped', StrictGenoCheck = TRUE, CalcLLR = TRUE, Tfilter = -10, Tassign = 0)
gmr10 <- GetMaybeRel(GenoM = filtered_gmmat_po_check_10)
# 8 trios. valid with hiphop.
################################################################################


## trying custom functions
################################################################################

# i have now customized the GetMaybeRel function using gmr_custom.R
# now, you can pass in the MaxMismatch argument to relax the ME and OH filters to be the same number.
# this may need to be further refined, but let's try it...

# ---- call custom function with relaxed OH/ME ---- #
gmr_custom <- GetMaybeRel_Custom(
  GenoM = filtered_gmmat_po_check,
  Module = "par",     
  MaxMismatch = 100,   # <- RELAX threshold here
  quiet = FALSE,
  Tassign = 0, 
  Tfilter = -10,
  MaxPairs = 7*nrow(filtered_gmmat_po_check)
)
#48 po pairs, 1 po trio

# retain 50% of the original loci
samp_loci <- sample(ncol(filtered_gmmat_po_check), 648)
filtered_gmmat_po_check_50 <- filtered_gmmat_po_check[, samp_loci]
gmr_custom_50 <- GetMaybeRel_Custom(
  GenoM = filtered_gmmat_po_check_50,
  Module = "par",     
  MaxMismatch = 648,   # <- RELAX threshold here
  quiet = FALSE,
  Tassign = 0, 
  Tfilter = -10,
  MaxPairs = 7*nrow(filtered_gmmat_po_check)
)
# got only 46 po pairs that time...


#let's try 10%
samp_loci <- sample(ncol(filtered_gmmat_po_check), 129)
filtered_gmmat_po_check_10 <- filtered_gmmat_po_check[, samp_loci]
outfull <- sequoia_custom(GenoM = filtered_gmmat_po_check_10, 
                          Module = 'ped',
                          MaxMismatch = 129,
                          Tassign = 0.01,
                          Tfilter = -100,
                          StrictGenoCheck = TRUE, 
                          CalcLLR = TRUE)

gmr_custom_10 <- GetMaybeRel_Custom(
  GenoM = filtered_gmmat_po_check_10,
  Module = "par",
  MaxMismatch = 129,   # <- RELAX threshold here
  quiet = FALSE,
  Tassign = 0.01, 
  Tfilter = -100,
  MaxPairs = 7*nrow(filtered_gmmat_po_check)
)
# same w this one.
# even including the custom functions!!!
################################################################################



## does including LH data help? 
################################################################################
setwd("/Users/samjohnson/Documents/Sauger_042225/GeneticData/Sequoia/Sequoia_Inp/")
LH_Data <- read.csv(file = "testindivs_LH.csv", header = TRUE)
# alright, what's going on here with the missing individual?
missingind <- testsamp$sample[!testsamp$sample %in% LH_Data$ID]
# first thing is to check the extractions, readme, and hiphop scripts.
# i understand 6757 wasn't sequenced, but what's up with 6436. it's on the plate map!
# 6436 wasn't spawned. 

# here's what we need to do:
# load in LH data
# change the F's to 1 and M's to 2
LH_Data$Sex[LH_Data$Sex == "M"] <- 2
LH_Data$Sex[LH_Data$Sex == "F"] <- 1
LH_Data$Sex <- as.numeric(LH_Data$Sex)
str(LH_Data)

LH_Data <- data.frame(LH_Data, BY.min = NA, BY.max = NA)
LH_Data$BY.min[1:95] <- 2005
LH_Data$BY.max[1:95] <- 2005
LH_Data$BY.min[96:nrow(LH_Data)] <- 2000
LH_Data$BY.max[96:nrow(LH_Data)] <- 2000

# filter the genotype matrix so that it only includes the ids from LH_Data$ID
filtered_gmmat_po_check <- filtered_gmmat_po_check[rownames(filtered_gmmat_po_check) %in% LH_Data$ID, , drop = FALSE]
dim(filtered_gmmat_po_check) # 208 indivs, 1296 loci

# run REGULAR SEQUOIA on that with the LH_Data specified.
# LifeHistData = LH_Data (make sure it's a df)
outfull <- sequoia(GenoM = filtered_gmmat_po_check_10, 
                   Module = 'par',
                   LifeHistData = LH_Data,
                   MaxMismatch = 129,
                   Tassign = 0.01,
                   Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)

gmr_reg_LHinc <- GetMaybeRel(GenoM = filtered_gmmat_po_check,
                             Module = "par",
                             MaxMismatch = 129,   # <- RELAX threshold here
                             LifeHistData = LH_Data,
                             quiet = FALSE,
                             Tassign = 0.01, 
                             Tfilter = -100,
                             MaxPairs = 7*nrow(filtered_gmmat_po_check))

# let's try with 75% loci
samp_loci <- sample(ncol(filtered_gmmat_po_check), 972)
filtered_gmmat_po_check_75 <- filtered_gmmat_po_check[ ,samp_loci]
outfull_reg_75 <- sequoia(GenoM = filtered_gmmat_po_check_75, 
                   Module = 'ped',
                   LifeHistData = LH_Data,
                   MaxMismatch = 129,
                   Tassign = 0.01,
                   Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)

outfull_custom_75 <- sequoia_custom(GenoM = filtered_gmmat_po_check_75, 
                                    Module = 'ped',
                                    LifeHistData = LH_Data,
                                    MaxMismatch = 129,
                                    Tassign = 0.01,
                                    Tfilter = -100,
                                    StrictGenoCheck = TRUE, 
                                    CalcLLR = TRUE)

gmr_reg_LHinc <- GetMaybeRel(GenoM = filtered_gmmat_po_check_75,
                             #SeqList = outfull_75,
                             Module = "ped",
                             MaxMismatch = 129,   # <- RELAX threshold here
                             LifeHistData = LH_Data,
                             quiet = FALSE,
                             Tassign = 0.01, 
                             Tfilter = -100,
                             MaxPairs = 7*nrow(filtered_gmmat_po_check))

gmr_custom_LHinc <- GetMaybeRel_Custom(GenoM = filtered_gmmat_po_check_75,
                                       #SeqList = outfull_75,
                                       Module = "ped",
                                       MaxMismatch = 129,   # <- RELAX threshold here
                                       LifeHistData = LH_Data,
                                       quiet = FALSE,
                                       Tassign = 0.01, 
                                       Tfilter = -100,
                                       MaxPairs = 7*nrow(filtered_gmmat_po_check))

# let's try with 40% loci
  # 75% was a reasonable test, but let's get into a good number of loci
  # and start messing with combinations of ped/par, including/excluding pedigree for gmr,
  # and reg/custom functions
samp_loci <- sample(ncol(filtered_gmmat_po_check), 518)
filtered_gmmat_po_check_40 <- filtered_gmmat_po_check[ ,samp_loci]
outfull_reg_40 <- sequoia(GenoM = filtered_gmmat_po_check_40, 
                          Module = 'par',
                          LifeHistData = LH_Data,
                          MaxMismatch = 129,
                          Tassign = 0.01,
                          Tfilter = -100,
                          StrictGenoCheck = TRUE, 
                          CalcLLR = TRUE)
# PED: assigned 89 dams and 87 sires to 208 + 47 individuals (real + dummy)
# PAR: assigned 25 dams and 22 sires to 208 individuals

outfull_custom_40 <- sequoia_custom(GenoM = filtered_gmmat_po_check_40, 
                                    Module = 'par',
                                    LifeHistData = LH_Data,
                                    MaxMismatch = 129,
                                    Tassign = 0.01,
                                    Tfilter = -100,
                                    StrictGenoCheck = TRUE, 
                                    CalcLLR = TRUE)
# PED: assigned 83 dams and 81 sires to 208 + 51 individuals (real + dummy)
# PAR: assigned 20 dams and 11 sires to 208 individuals

gmr_reg_LHinc <- GetMaybeRel(GenoM = filtered_gmmat_po_check_75,
                             SeqList = outfull_reg_40,
                             Module = "par",
                             MaxMismatch = 129,   # <- RELAX threshold here
                             LifeHistData = LH_Data,
                             quiet = FALSE,
                             Tassign = 0.01, 
                             Tfilter = -100,
                             MaxPairs = 7*nrow(filtered_gmmat_po_check))
# PED/condition on outfull_reg_40: Found 9 likely parent-offspring pairs, and 335, other non-assigned pairs of possible relatives
# PED/NOT CONDITIONING:  Found 42 likely parent-offspring pairs, and 382, other non-assigned pairs of possible relatives

# PAR/condition on outfull_reg_40: Found 4 likely parent-offspring pairs, and 1, other non-assigned pairs of possible relatives
# PAR/NOT CONDITIONING: Found 34 likely parent-offspring pairs, and 4, other non-assigned pairs of possible relatives


gmr_custom_LHinc <- GetMaybeRel_Custom(GenoM = filtered_gmmat_po_check_75,
                                       SeqList = outfull_75,
                                       Module = "par",
                                       MaxMismatch = 129,   # <- RELAX threshold here
                                       LifeHistData = LH_Data,
                                       quiet = FALSE,
                                       Tassign = 0.01, 
                                       Tfilter = -100,
                                       MaxPairs = 7*nrow(filtered_gmmat_po_check))
# PED/condition on outfull_reg_40: 4 likely parent-offspring pairs, and 1, other non-assigned pairs of possible relatives
# PED/NOT CONDITIONING: Found 42 likely parent-offspring pairs, and 382 other non-assigned pairs of possible relatives
                        # Found 1 parent-parent-offspring trios

# PAR/condition on outfull_reg_40: Found 4 likely parent-offspring pairs, and 3 other non-assigned pairs of possible relatives
# PAR/NOT CONDITIONING: Found 46 likely parent-offspring pairs, and 14 other non-assigned pairs of possible relatives
                        # Found 1 parent-parent-offspring trios
################################################################################




## heard from Nancy Chen that sequoia has trouble with ANY error in large datasets
## does a min mean depth filter of 15 tighten it up?
################################################################################
setwd("/Users/samjohnson/Documents/Sauger_042225/GeneticData/Sequoia/Sequoia_Inp/WAE_bwamem_hardfilter")

# read in genotype matrix
mat <- read.table(file = "hard_variants_pflav_mem_t2_bial_noindels_q20_mindep15_maxdep75_maf30_miss95.012_conv", header = FALSE, sep = "\t", 
                  na.strings = c("NA", "-1"))

# correct row names as sample id's
mat<- mat[, -1]
gmmat <- as.matrix(mat)
ind <- read.table(file = "hard_variants_pflav_mem_t2_bial_noindels_q20_mindep15_maxdep75_maf30_miss95.012.indv", header = FALSE)

str(ind)
ind <- ind %>% 
  rename(sample = V1)
rownames(gmmat) <- ind$sample

# read in parent offspring (F0 and Test F1 sample id's)
testsamp <- read.csv(file = "posampleids.csv", header = TRUE) # 210 samples, 114 Parents, 96 Test F1's

# read in lh data
LH_Data <- read.csv(file = "testindivs_LH.csv", header = TRUE)
# alright, what's going on here with the missing individual?
missingind <- testsamp$sample[!testsamp$sample %in% LH_Data$ID]
# first thing is to check the extractions, readme, and hiphop scripts.
# i understand 6757 wasn't sequenced, but what's up with 6436. it's on the plate map!
# 6436 wasn't spawned. 

# here's what we need to do:
# change the F's to 1 and M's to 2
LH_Data$Sex[LH_Data$Sex == "M"] <- 2
LH_Data$Sex[LH_Data$Sex == "F"] <- 1
LH_Data$Sex <- as.numeric(LH_Data$Sex)
str(LH_Data)

# add birth year info
LH_Data <- data.frame(LH_Data, BY.min = NA, BY.max = NA)
LH_Data$BY.min[1:95] <- 2005
LH_Data$BY.max[1:95] <- 2005
LH_Data$BY.min[96:nrow(LH_Data)] <- 2000
LH_Data$BY.max[96:nrow(LH_Data)] <- 2000
# NEED TO WRITE THIS OUT FOR FUTURE USE

# filter the genotype matrix so that it only includes the ids from LH_Data$ID
gmmat_po <- gmmat[rownames(gmmat) %in% LH_Data$ID, , drop = FALSE]
dim(gmmat) # 208 indivs, 978 loci

# check for all heterozygous sites since those likely will not be informative.
all_ones <- apply(gmmat_po_check, 2, function(col) all(col == 1))
filtered_gmmat_po_check <- gmmat_po_check[, !all_ones]
# removed two sites

# check genotype matrix for samples/loci to be excluded
filtered_gmmat_po_check <- CheckGeno(filtered_gmmat_po_check, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))    
# 95 Test F1's, 114 Parents and 976 SNPs going into this.

# run sequoia on the checked genotype matrix
outfull <- sequoia_custom(GenoM = filtered_gmmat_po_check, 
                   Module = 'par',
                   LifeHistData = LH_Data,
                   MaxMismatch = 978,
                   Tassign = 0.01,
                   Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)
# adding all of this stuff makes it WORSE GOD DAMMIT

gmr<- GetMaybeRel_Custom(GenoM = filtered_gmmat_po_check_75,
                             SeqList = outfull,
                             Module = "par",
                             MaxMismatch = 796,   
                             LifeHistData = LH_Data,
                             quiet = FALSE,
                             Tassign = 0.01, 
                             Tfilter = -100,
                             MaxPairs = 7*nrow(filtered_gmmat_po_check))
# alright there seems to be no rhyme or reason to how many i can add and i'm starting
# to freak out so i'm going to call it here and wait to hear from jisca.














