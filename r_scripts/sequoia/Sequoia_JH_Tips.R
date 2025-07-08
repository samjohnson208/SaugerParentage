## Sequoia_JH_Tips.R by SPJ 070325
## PURPOSE: to understand problems assoc. with Sequoia after getting a wonderful
# list of suggestions from Dr. Jisca Huisman. In order of implementation, they are:
# (numbered by the scheme from her email)
# 3. filter the vcfs for LD by thinnning and retry.
# 1. Increase the Tassign threshold
# 2. Alter the complex argument from sequoia()
# 4. Try an older sequoia program version (currently running 2.11.2)
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

setwd("/Users/samjohnson/Documents/Sauger_042225/GeneticData/Sequoia/Sequoia_Inp/svit_mem_thinned")

# read in genotype matrix
mat <- read.table(file = "hard_variants_svit_mem_bial_noindels_q20_mindep15_maxdep75_maf30_miss95_thin10K.012_conv", 
                  header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

# correct row names as sample id's
mat <- mat[, -1]
gmmat <- as.matrix(mat)
ind <- read.table(file = "hard_variants_svit_mem_bial_noindels_q20_mindep15_maxdep75_maf30_miss95_thin10K.012.indv", header = FALSE)

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
dim(gmmat_po) # 208 indivs, 865 loci

# check for all heterozygous sites since those likely will not be informative.
all_ones <- apply(gmmat_po, 2, function(col) all(col == 1))
gmmat_po_check <- gmmat_po[, !all_ones]
dim(gmmat_po_check)
# lost one snp

# check genotype matrix for samples/loci to be excluded
filtered_gmmat_po_check <- CheckGeno(gmmat_po_check, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))    
# 95 Test F1's, 114 Parents and 864 SNPs. 

# run sequoia on the checked genotype matrix
outfull <- sequoia(GenoM = filtered_gmmat_po_check, 
                   Module = 'par',
                   LifeHistData = LH_Data,
                 # MaxMismatch = 978,
                   Tassign = 1.0,
                 # Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)

gmr<- GetMaybeRel(GenoM = filtered_gmmat_po_check,
                  SeqList = outfull,
                  Module = "par",
                # MaxMismatch = 796,   
                  LifeHistData = LH_Data,
                  quiet = FALSE,
                  Tassign = 1.0, 
                # Tfilter = -100,
                  MaxPairs = 7*nrow(filtered_gmmat_po_check))

# 9 pairs. LH Data included. Too many loci.
# Sample some loci, try the same. Then Tassign. Then Complex. Then older version.

# Increasing Tassign to 1.0 got me to 17 with 864 SNPs. Going to start dropping to 
# some reduced numbers of SNPs.

## ~~~~~~~~~~ 600 SNPs ~~~~~~~~~~ ##

samp_600 <- sample(ncol(filtered_gmmat_po_check), 600)
filtered_gmmat_po_check_600 <- filtered_gmmat_po_check[, samp_600]

# run sequoia on the checked genotype matrix
outfull <- sequoia(GenoM = filtered_gmmat_po_check_600, 
                   Module = 'par',
                   LifeHistData = LH_Data,
                   # MaxMismatch = 978,
                   Tassign = 1.0,
                   # Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)

gmr<- GetMaybeRel(GenoM = filtered_gmmat_po_check_600,
                  #SeqList = outfull,
                  Module = "par",
                  # MaxMismatch = 796,   
                  LifeHistData = LH_Data,
                  quiet = FALSE,
                  Tassign = 1.0, 
                  # Tfilter = -100,
                  MaxPairs = 7*nrow(filtered_gmmat_po_check))
# 18 PO pairs, 6 other non-assigned pairs of possible relatives.
# 24 PO pairs, 0 other non-assigned pairs of possible relatives.
# Likely pretty dependent on which loci are selected. Keep reducing.
# When you don't condition on the pedigree from outfull you do WAY better.
# 92 likely PO pairs, 10 others.

## ~~~~~~~~~~ 500 SNPs ~~~~~~~~~~ ##

samp_500 <- sample(ncol(filtered_gmmat_po_check), 500)
filtered_gmmat_po_check_500 <- filtered_gmmat_po_check[, samp_500]

# run sequoia on the checked genotype matrix
outfull <- sequoia(GenoM = filtered_gmmat_po_check_500, 
                   Module = 'par',
                   LifeHistData = LH_Data,
                   # MaxMismatch = 978,
                   Tassign = 1.0,
                   # Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)

gmr <- GetMaybeRel(GenoM = filtered_gmmat_po_check_500,
                  # SeqList = outfull,
                  Module = "ped",
                  # MaxMismatch = 796,   
                  LifeHistData = LH_Data,
                  quiet = FALSE,
                  Tassign = 1.0, 
                  # Tfilter = -100,
                  MaxPairs = 7*nrow(filtered_gmmat_po_check))
# 17 PO pairs, 5 other non-assigned pairs of possible relatives. (Module = "par")
# 17 PO pairs, 5 other non-assigned pairs of possible relatives. (Module = "par")
# When you don't condition on the pedigree from outfull you do WAY better.
# 93 likely PO pairs, 10 others.
# Resample, 116 PO pairs, 5 others.
# What if you don't specify the module?
# 116 PO pairs, 5 others, 9 trios with module commented out. 
# (I guess that means that default module is "par", but you get trios if you don't specify?)
# 112 PO pairs and 319 other. 9 trios with Module = "ped".
# I think this is the largest amount of assignments I may have ever gotten.


## ~~~~~~~~~~ 400 SNPs ~~~~~~~~~~ ##

samp_400 <- sample(ncol(filtered_gmmat_po_check), 400)
filtered_gmmat_po_check_400 <- filtered_gmmat_po_check[, samp_400]

# run sequoia on the checked genotype matrix
outfull <- sequoia(GenoM = filtered_gmmat_po_check_400, 
                   Module = 'par',
                   LifeHistData = LH_Data,
                   # MaxMismatch = 978,
                   Tassign = 1.0,
                   # Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)

gmr <- GetMaybeRel(GenoM = filtered_gmmat_po_check_400,
                   # SeqList = outfull,
                   Module = "ped",
                   # MaxMismatch = 796,   
                   LifeHistData = LH_Data,
                   quiet = FALSE,
                   Tassign = 1.0, 
                   # Tfilter = -100,
                   MaxPairs = 7*nrow(filtered_gmmat_po_check))
# with SeqList specified and module = "par": 25 PO pairs, 0 others, 2 trios.
# seqlist commented out: 110 pairs, 11 others, 10 trios.
# module = "ped": 114 PO pairs, 335 others, 10 trios.
# minimal differences with resampling at this point

## ~~~~~~~~~~ 300 SNPs ~~~~~~~~~~ ##

samp_300 <- sample(ncol(filtered_gmmat_po_check), 300)
filtered_gmmat_po_check_300 <- filtered_gmmat_po_check[, samp_300]

# run sequoia on the checked genotype matrix
outfull <- sequoia(GenoM = filtered_gmmat_po_check_300, 
                   # Module = 'par',
                   LifeHistData = LH_Data,
                   # MaxMismatch = 978,
                   Tassign = 1.0,
                   # Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)

gmr <- GetMaybeRel(GenoM = filtered_gmmat_po_check_300,
                   # SeqList = outfull,
                   Module = "ped",
                   # MaxMismatch = 796,   
                   LifeHistData = LH_Data,
                   quiet = FALSE,
                   Tassign = 1.0,
                   # Tfilter = -100
                   MaxPairs = 7*nrow(filtered_gmmat_po_check))
# with SeqList specified and module = "par": 20 PO pairs, 126 others, no trios.
# seqlist commented out: 98 pairs pairs, 22 others, 17 trios.
# module = "ped": 94 pairs pairs, 344 others others, 17 trios.
# resampling, module par: 115 PO pairs, 14 others, 21 trios. 
# resampling, module ped: 106 PO pairs, 320 others, 20 trios.

## ~~~~~~~~~~ 200 SNPs ~~~~~~~~~~ ##

samp_200 <- sample(ncol(filtered_gmmat_po_check), 200)
filtered_gmmat_po_check_200 <- filtered_gmmat_po_check[, samp_200]

# run sequoia on the checked genotype matrix
outfull <- sequoia(GenoM = filtered_gmmat_po_check_200, 
                   # Module = 'par',
                   LifeHistData = LH_Data,
                   # MaxMismatch = 978,
                   Tassign = 1.0,
                   # Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)

gmr <- GetMaybeRel(GenoM = filtered_gmmat_po_check_200,
                   # SeqList = outfull,
                   Module = "ped",
                   # MaxMismatch = 796,   
                   LifeHistData = LH_Data,
                   quiet = FALSE,
                   Tassign = 1.0,
                   # Tfilter = -100
                   MaxPairs = 7*nrow(filtered_gmmat_po_check))
# par, SeqList commented out: 84 pairs, 20 others, 12 trios.
# ped, SeqList commented out: 73 pairs, 336 others, 10 trios.
# resampling, ped, SeqList commented out: 95 pairs, 345 others, 25 trios.

## ~~~~~~~~~~ 100 SNPs ~~~~~~~~~~ ##

samp_100 <- sample(ncol(filtered_gmmat_po_check), 100)
filtered_gmmat_po_check_100 <- filtered_gmmat_po_check[, samp_100]

# run sequoia on the checked genotype matrix
outfull <- sequoia(GenoM = filtered_gmmat_po_check_100, 
                   # Module = 'par',
                   LifeHistData = LH_Data,
                   # MaxMismatch = 978,
                   Tassign = 1.0,
                   # Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)

gmr <- GetMaybeRel(GenoM = filtered_gmmat_po_check_100,
                   # SeqList = outfull,
                   Module = "ped",
                   # MaxMismatch = 796,   
                   LifeHistData = LH_Data,
                   quiet = FALSE,
                   Tassign = 1.0,
                   # Tfilter = -100
                   MaxPairs = 7*nrow(filtered_gmmat_po_check))
# par, SeqList: 51 pairs, 11 others, 9 trios.
# par, SeqList commented out: 117 pairs, 24 others, 33 trios.
# ped, SeqList commented out: 104 pairs, 234 others, 33 trios.
# resampling, ped, SeqList commented out: 102 pairs, 284 others, 39 trios.
# most pairs here for sure with 39. resampling causes pairs to bounce around in the
# 20's and 30's. makes me wonder what else is contributing to the quality of these snps.

# Alright since we kept the most trios here, and just about the most pairs, I'm going
# to go ahead and introduce the Complex argument here. 

samp_100 <- sample(ncol(filtered_gmmat_po_check), 100)
filtered_gmmat_po_check_100 <- filtered_gmmat_po_check[, samp_100]

# run sequoia on the checked genotype matrix
outfull <- sequoia(GenoM = filtered_gmmat_po_check_100, 
                   # Module = 'par',
                   LifeHistData = LH_Data,
                   # MaxMismatch = 978,
                   Complex = "full",
                   Tassign = 1.0,
                   # Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)

gmr <- GetMaybeRel(GenoM = filtered_gmmat_po_check_100,
                   # SeqList = outfull,
                   Module = "ped",
                   # MaxMismatch = 796,
                   Complex = "simp",
                   LifeHistData = LH_Data,
                   quiet = FALSE,
                   Tassign = 1.0,
                   # Tfilter = -100
                   MaxPairs = 7*nrow(filtered_gmmat_po_check_100))
# ped, SeqList = outfull, Complex = "full": 38 pairs, 93 others, 3 trios
# ped, SeqList commented out, Complex = "full": 97 pairs, 248 others, 21 trios
# par, SeqList commented out, Complex = "full": 120 pairs, 25 others, 25 trios
# resamp: par, SeqList commented out, Complex = "full": 141 pairs, 17 others, 40 trios
# resamp: par, SeqList commented out, Complex = "full": 510 pairs, 11 others, 41 trios
# par, SeqList commented out, Complex = "simp": 150 pairs, 19 others, 45 trios. WOW.
# par, SeqList = outfull, Complex = "simp": 63, 0 others, 12 trios
# ped, SeqList = outfull, Complex = "simp": 56 pairs, 121 others, 12 trios.
# ped, SeqList commented out, Complex = "simp": 

gmr <- GetMaybeRel_Custom(GenoM = filtered_gmmat_po_check_100,
                        # SeqList = outfull,
                          Module = "par",
                          MaxMismatch = 100,
                          Complex = "simp",
                          LifeHistData = LH_Data,
                          quiet = FALSE,
                          Tassign = 1.0,
                        # Tfilter = -100
                          MaxPairs = 7*nrow(filtered_gmmat_po_check))

# custom: 150 pairs, 11 others, 25 trios. hmmmmm.... Still ZERO mismatches. 
# Try the addition of more mismatches? Tried with MaxMismatch = 100 and there were still
# no trio assignments where the parent had more than one OH mismatch, or more than 1 ME
# mismatch.

################################################################################
#### RUNNING ON FURTHER THINNED (LD FILTERED) FILES
################################################################################

setwd("/Users/samjohnson/Documents/Sauger_042225/GeneticData/Sequoia/Sequoia_Inp/svit_mem_thinned")

# read in genotype matrix
mat1M <- read.table(file = "hard_variants_svit_mem_bial_noindels_q20_mindep15_maxdep75_maf30_miss95_thin1M.012_conv", 
                  header = FALSE, sep = "\t", na.strings = c("NA", "-1"))
mat2.5M <- read.table(file = "hard_variants_svit_mem_bial_noindels_q20_mindep15_maxdep75_maf30_miss95_thin2.5M.012_conv", 
                    header = FALSE, sep = "\t", na.strings = c("NA", "-1"))
mat5M <- read.table(file = "hard_variants_svit_mem_bial_noindels_q20_mindep15_maxdep75_maf30_miss95_thin5M.012_conv", 
                    header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

# correct row names as sample id's
mat1M <- mat1M[, -1]
mat2.5M <- mat2.5M[, -1]
mat5M <- mat5M[, -1]

mat1M <- as.matrix(mat1M)
mat2.5M <- as.matrix(mat2.5M)
mat5M <- as.matrix(mat5M)

ind1M <- read.table(file = "hard_variants_svit_mem_bial_noindels_q20_mindep15_maxdep75_maf30_miss95_thin1M.012.indv", header = FALSE)
ind2.5M <- read.table(file = "hard_variants_svit_mem_bial_noindels_q20_mindep15_maxdep75_maf30_miss95_thin2.5M.012.indv", header = FALSE)
ind5M <- read.table(file = "hard_variants_svit_mem_bial_noindels_q20_mindep15_maxdep75_maf30_miss95_thin5M.012.indv", header = FALSE)

ind1M <- ind1M %>% 
  rename(sample = V1)
ind2.5M <- ind2.5M %>% 
  rename(sample = V1)
ind5M <- ind5M %>% 
  rename(sample = V1)

rownames(mat1M) <- ind1M$sample
rownames(mat2.5M) <- ind2.5M$sample
rownames(mat5M) <- ind5M$sample

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
mat1M <- mat1M[rownames(mat1M) %in% LH_Data$ID, , drop = FALSE]
dim(mat1M) # 208 indivs, 371 loci

mat2.5M <- mat2.5M[rownames(mat2.5M) %in% LH_Data$ID, , drop = FALSE]
dim(mat2.5M) # 208 indivs, 220 loci

mat5M <- mat5M[rownames(mat5M) %in% LH_Data$ID, , drop = FALSE]
dim(mat5M) # 208 indivs, 140 loci

# check for all heterozygous sites since those likely will not be informative.
all_ones_1M <- apply(mat1M, 2, function(col) all(col == 1))
all_ones_2.5M <- apply(mat2.5M, 2, function(col) all(col == 1))
all_ones_5M <- apply(mat5M, 2, function(col) all(col == 1))

mat1M <- mat1M[, !all_ones_1M]
mat2.5M <- mat2.5M[, !all_ones_2.5M]
mat5M <- mat5M[, !all_ones_5M]

dim(mat1M)
# lost one snp
dim(mat2.5M)
# kept all
dim(mat5M)
# kept all

# check genotype matrix for samples/loci to be excluded
check1M <- CheckGeno(mat1M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))    
# 95 Test F1's, 113 Parents and 370 SNPs

check2.5M <- CheckGeno(mat2.5M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))    
# 95 Test F1's, 113 Parents and 220 SNPs

check5M <- CheckGeno(mat5M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))    
# 95 Test F1's, 113 Parents and 140 SNPs

## ~~~~~~~~~~ 1Mbp, 370 SNPs ~~~~~~~~~~ ##

# run sequoia on the checked genotype matrix
outfull <- sequoia(GenoM = check1M, 
                   Module = 'par',
                   LifeHistData = LH_Data,
                   # MaxMismatch = 978,
                   Complex = "simp",
                   Tassign = 1.0,
                   # Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)
# par, complex="full": 78 dams, 71 sires to 208 indivs...
# par, complex="simp": 81 dams 77 sires to 208 indivs...

gmr <- GetMaybeRel(GenoM = check1M,
                   # SeqList = outfull,
                   Module = "par",
                   # MaxMismatch = 796,
                   Complex = "simp",
                   LifeHistData = LH_Data,
                   quiet = FALSE,
                   Tassign = 1.0,
                   # Tfilter = -100
                   MaxPairs = 7*nrow(check1M))
# par, complex="simp", SeqList = outfull: 8 pairs, 0 others, 0 trios.
# par, complex="simp", no SeqList: 120 pairs, 29 others, 15 trios.

## ~~~~~~~~~~ 2.5Mbp, 220 SNPs ~~~~~~~~~~ ##

outfull <- sequoia(GenoM = check2.5M, 
                   Module = 'par',
                   LifeHistData = LH_Data,
                   # MaxMismatch = 978,
                   Complex = "simp",
                   Tassign = 1.0,
                   # Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)
# par, simp: 45 dams, 42 sires. 208 indivs.

gmr <- GetMaybeRel(GenoM = check2.5M,
                   # SeqList = outfull,
                   Module = "par",
                   # MaxMismatch = 796,
                   Complex = "simp",
                   LifeHistData = LH_Data,
                   quiet = FALSE,
                   Tassign = 1.0,
                   # Tfilter = -100
                   MaxPairs = 7*nrow(check2.5M))
# par, seqlist = outfull, simp: 31 pairs, 0 others, 3 trios.
# par, no seqlist, simp: 101 pairs, 54 others, 26 trios.

## ~~~~~~~~~~ 5Mbp, 140 SNPs ~~~~~~~~~~ ##

outfull <- sequoia(GenoM = check5M, 
                   Module = 'par',
                   LifeHistData = LH_Data,
                   # MaxMismatch = 978,
                   Complex = "simp",
                   Tassign = 1.0,
                   # Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)
# par, simp: 36 dams 40 sires 208 indivs. 

gmr <- GetMaybeRel(GenoM = check5M,
                   # SeqList = outfull,
                   Module = "par",
                   # MaxMismatch = 796,
                   Complex = "simp",
                   LifeHistData = LH_Data,
                   quiet = FALSE,
                   Tassign = 1.0,
                   # Tfilter = -100,
                   MaxPairs = 7*nrow(check5M))
# par, SeqList = outfull, simp: 47 po, 0 others, 6 trios.
# par, no seqlist, simp: 108 po, 47 others, 27 trios. 

# take 100 of the 5M snps
samp_100 <- sample(ncol(check5M), 100)
check5M_100 <- check5M[, samp_100]

gmr <- GetMaybeRel(GenoM = check5M_100,
                   # SeqList = outfull,
                   Module = "par",
                   # MaxMismatch = 796,
                   Complex = "simp",
                   LifeHistData = LH_Data,
                   quiet = FALSE,
                   Tassign = 1.0,
                   # Tfilter = -100,
                   MaxPairs = 7*nrow(check5M_100))
# 131 pairs, 36 others, 34 trios... hmmm... Wonder what's up with these...
# Obviously a discrepancy in the quality of some of these SNPs. Wonder if it has
# to do with the allele frequencies. You'd think that a mix of SNPs with varying
# MAFs would be good but I suppose not if the recommended MAF filter is 30% and up.

# try with error?
error_rate <- c(0.0001, 0.0001, 0.0001)
gmr <- GetMaybeRel(GenoM = check5M_100,
                   Err = error_rate,
                   # SeqList = outfull,
                   Module = "par",
                   # MaxMismatch = 796,
                   Complex = "simp",
                   LifeHistData = LH_Data,
                   quiet = FALSE,
                   Tassign = 1.0,
                   # Tfilter = -100,
                   MaxPairs = 7*nrow(check5M_100))
# WITH ERROR 0.01 WE GET 164 PAIRS, 49 OTHERS, 80 TRIOS.
# WITH ERROR 0.02 WE GET 170 PAIRS, 64 OTHERS, 86 TRIOS.
# WITH ERROR 0.025 WE GET 171 PAIRS, 73 OTHERS, 91 TRIOS.
# WITH ERROR RATE 0.0001 132 PAIRS, 36 OTHERS, 35 TRIOS.

# additional things to try: filter for sites with NO missing data
# see how many we pull, see if increasing the number of sites, but having a really
# small error rate changes things. lots of dials to turn.


