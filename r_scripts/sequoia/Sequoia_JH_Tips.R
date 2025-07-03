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
                 # Tassign = 0.01,
                 # Tfilter = -100,
                   StrictGenoCheck = TRUE, 
                   CalcLLR = TRUE)

gmr<- GetMaybeRel(GenoM = filtered_gmmat_po_check,
                  SeqList = outfull,
                  Module = "par",
                # MaxMismatch = 796,   
                  LifeHistData = LH_Data,
                  quiet = FALSE,
                # Tassign = 0.01, 
                # Tfilter = -100,
                  MaxPairs = 7*nrow(filtered_gmmat_po_check))
# 9 pairs. LH Data included. Too many loci.
# Sample some loci, try the same. Then Tassign. Then Complex. Then older version.









