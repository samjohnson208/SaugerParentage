## Sequoia_ContamFilt_mindep8_md5.R by SPJ 091225
## PURPOSE: I have compiled a data of sequoia results for various n snps, md, and
# error. (sequoia_params.xlsx). all this comes from a dataset that has been
# cleaned using illumina, phix, and ecoli databases, and has been cleaned using
# FASTP (remove bases with q<20, reads shorter than 50bp, poly g tails, and any 
# remaining sequence adapters). This dataset was then aligned to the walleye reference
# using BWA MEM, variant calling was run, vcf was reheadered, and the data were 
# filtered using the slurm_1 and slurm_3 filtering scripts for the following params.

# rawfiltered:see bcftools mpileup in slurm_sauger_variants.sh
# contam - contaminant filtering mentioned above
# fastp - see above
# bial - keep only biallelic sites
# no indels - keep sites that are not insertions or deletions
# q40 - keep sites with site quality > 40
# mindep8 - minimum mean read depth 4 or above (removes low quality sites, low confidence calls)
# maxdep75 - maximum mean read depth 75 or lower (removes paralogs)
# maf30 - minor allele frequency 30% or higher
# miss95 - each snp needs to have 5% or less missing data

# and thinned to various degrees to filter for LD (25K, 50K, 75K, 100K to 900K, 
# we already ran miss 95 1M to 5M in Sequoia_Contam_mindep8.R).

# these vcf -> .012 came from me wanting to be more confident in my genotype calls
# which is why site qual went up from 20 to 40, and mindep went up to 8. 

# we already did lower proportions of missing data with 40 - 400+ snps, and we
# saw that more sites, with a bit more md, and a bit more error, gave me the highest
# assignment and accuracy rates for the test f1's.

# this is me working to CONTINUE to fill the sequoia_params table. we want a spread
# of sites that's the same with various degrees of md and thinning. goal is to expand
# the 5% and 10% md portions to include a wide range of thinning and a wide range
# of nsnps in order to have valid comparisons and to explore the range that sequoia
# has across these axes.

## libraries
################################################################################

# install.packages("sequoia")
# 65
library(sequoia)

# install.packages("dplyr")
# 65
library(dplyr)

################################################################################

setwd("/Users/samjohnson/Documents/Sauger_042225/GeneticData/Sequoia/Sequoia_Inp/contam_fastp_svit_mem/firstfilt_hardfilt_thinned/mindep8_maf30/geno_mat")

################################################################################
#### RUNNING ON FURTHER THINNED (LD FILTERED) FILES ####
################################################################################

##### ---- Prepping Matrices ---- #####

mat_thin25K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin25K.012_conv",
                          header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin50K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin50K.012_conv",
                          header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin75K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin75K.012_conv",
                          header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin100K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012_conv",
                           header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin200K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin200K.012_conv",
                           header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin300K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin300K.012_conv",
                           header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin400K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin400K.012_conv",
                           header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin500K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin500K.012_conv",
                           header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin600K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin600K.012_conv",
                           header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin700K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin700K.012_conv",
                           header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin800K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin800K.012_conv",
                           header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin900K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin900K.012_conv",
                           header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

# remove first column. ascending numbers 1 - 1184
mat_thin25K <- mat_thin25K[, -1]
mat_thin50K <- mat_thin50K[, -1]
mat_thin75K <- mat_thin75K[, -1]
mat_thin100K <- mat_thin100K[, -1]
mat_thin200K <- mat_thin200K[, -1]
mat_thin300K <- mat_thin300K[, -1]
mat_thin400K <- mat_thin400K[, -1]
mat_thin500K <- mat_thin500K[, -1]
mat_thin600K <- mat_thin600K[, -1]
mat_thin700K <- mat_thin700K[, -1]
mat_thin800K <- mat_thin800K[, -1]
mat_thin900K <- mat_thin900K[, -1]

# turn them to matrices
mat_thin25K <- as.matrix(mat_thin25K)
mat_thin50K <- as.matrix(mat_thin50K)
mat_thin75K <- as.matrix(mat_thin75K)
mat_thin100K <- as.matrix(mat_thin100K)
mat_thin200K <- as.matrix(mat_thin200K)
mat_thin300K <- as.matrix(mat_thin300K)
mat_thin400K <- as.matrix(mat_thin400K)
mat_thin500K <- as.matrix(mat_thin500K)
mat_thin600K <- as.matrix(mat_thin600K)
mat_thin700K <- as.matrix(mat_thin700K)
mat_thin800K <- as.matrix(mat_thin800K)
mat_thin900K <- as.matrix(mat_thin900K)

# yes, i am aware that these are all the same. trying to be thorough. would hate 
# to mess up a matrix because i was lazy on a step like this.
ind25 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin25K.012.indv", header = FALSE)
ind50 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin50K.012.indv", header = FALSE)
ind75 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin75K.012.indv", header = FALSE)
ind1 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012.indv", header = FALSE)
ind2 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin200K.012.indv", header = FALSE)
ind3 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin300K.012.indv", header = FALSE)
ind4 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin400K.012.indv", header = FALSE)
ind5 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin500K.012.indv", header = FALSE)
ind6 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin600K.012.indv", header = FALSE)
ind7 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin700K.012.indv", header = FALSE)
ind8 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin800K.012.indv", header = FALSE)
ind9 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin900K.012.indv", header = FALSE)

# rename first column to something meaningful
ind25 <- ind25 %>% 
  rename(sample = V1)
ind50 <- ind50 %>% 
  rename(sample = V1)
ind75 <- ind75 %>% 
  rename(sample = V1)
ind1 <-  ind1 %>% 
  rename(sample = V1)
ind2 <-  ind2 %>% 
  rename(sample = V1)
ind3 <-  ind3 %>% 
  rename(sample = V1)
ind4 <-  ind4 %>% 
  rename(sample = V1)
ind5 <-  ind5 %>% 
  rename(sample = V1)
ind6 <-  ind6 %>% 
  rename(sample = V1)
ind7 <-  ind7 %>% 
  rename(sample = V1)
ind8 <-  ind8 %>% 
  rename(sample = V1)
ind9 <-  ind9 %>% 
  rename(sample = V1)

# establish sample identities in the geno_mat
rownames(mat_thin25K) <- ind25$sample
rownames(mat_thin50K) <- ind50$sample
rownames(mat_thin75K) <- ind75$sample
rownames(mat_thin100K) <- ind1$sample
rownames(mat_thin200K) <- ind2$sample
rownames(mat_thin300K) <- ind3$sample
rownames(mat_thin400K) <- ind4$sample
rownames(mat_thin500K) <- ind5$sample
rownames(mat_thin600K) <- ind6$sample
rownames(mat_thin700K) <- ind7$sample
rownames(mat_thin800K) <- ind8$sample
rownames(mat_thin900K) <- ind9$sample

# read in scaffolds and positions
pos25 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin25K.012.pos", header = FALSE)
pos50 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin50K.012.pos", header = FALSE)
pos75 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin75K.012.pos", header = FALSE)
pos1 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012.pos", header = FALSE)
pos2 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin200K.012.pos", header = FALSE)
pos3 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin300K.012.pos", header = FALSE)
pos4 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin400K.012.pos", header = FALSE)
pos5 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin500K.012.pos", header = FALSE)
pos6 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin600K.012.pos", header = FALSE)
pos7 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin700K.012.pos", header = FALSE)
pos8 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin800K.012.pos", header = FALSE)
pos9 <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin900K.012.pos", header = FALSE)

# create full positions by combining the two columns
pos25$position <- paste(pos25$V1, pos25$V2, sep = "_")
pos50$position <- paste(pos50$V1, pos50$V2, sep = "_")
pos75$position <- paste(pos75$V1, pos75$V2, sep = "_")
pos1$position <-  paste(pos1$V1, pos1$V2, sep = "_")
pos2$position <-  paste(pos2$V1, pos2$V2, sep = "_")
pos3$position <-  paste(pos3$V1, pos3$V2, sep = "_")
pos4$position <-  paste(pos4$V1, pos4$V2, sep = "_")
pos5$position <-  paste(pos5$V1, pos5$V2, sep = "_")
pos6$position <-  paste(pos6$V1, pos6$V2, sep = "_")
pos7$position <-  paste(pos7$V1, pos7$V2, sep = "_")
pos8$position <-  paste(pos8$V1, pos8$V2, sep = "_")
pos9$position <-  paste(pos9$V1, pos9$V2, sep = "_")

# set positions as column names of geno_mat
colnames(mat_thin25K) <- pos25$position
colnames(mat_thin50K) <- pos50$position
colnames(mat_thin75K) <- pos75$position
colnames(mat_thin100K) <- pos1$position
colnames(mat_thin200K) <- pos2$position
colnames(mat_thin300K) <- pos3$position
colnames(mat_thin400K) <- pos4$position
colnames(mat_thin500K) <- pos5$position
colnames(mat_thin600K) <- pos6$position
colnames(mat_thin700K) <- pos7$position
colnames(mat_thin800K) <- pos8$position
colnames(mat_thin900K) <- pos9$position









##### ---- Adding Additional Data, Filtering by Individuals ---- #####

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
table(LH_Data$Sex) # 53 female F0, 60 male F0, 95 unknown sex test F1

# filter the genotype matrices so that it only includes the ids from LH_Data$ID
mat_thin25K <- mat_thin25K[rownames(mat_thin25K) %in% LH_Data$ID, , drop = FALSE]
mat_thin50K <- mat_thin50K[rownames(mat_thin50K) %in% LH_Data$ID, , drop = FALSE]
mat_thin75K <- mat_thin75K[rownames(mat_thin75K) %in% LH_Data$ID, , drop = FALSE]
mat_thin100K <- mat_thin100K[rownames(mat_thin100K) %in% LH_Data$ID, , drop = FALSE]
mat_thin200K <- mat_thin200K[rownames(mat_thin200K) %in% LH_Data$ID, , drop = FALSE]
mat_thin300K <- mat_thin300K[rownames(mat_thin300K) %in% LH_Data$ID, , drop = FALSE]
mat_thin400K <- mat_thin400K[rownames(mat_thin400K) %in% LH_Data$ID, , drop = FALSE]
mat_thin500K <- mat_thin500K[rownames(mat_thin500K) %in% LH_Data$ID, , drop = FALSE]
mat_thin600K <- mat_thin600K[rownames(mat_thin600K) %in% LH_Data$ID, , drop = FALSE]
mat_thin700K <- mat_thin700K[rownames(mat_thin700K) %in% LH_Data$ID, , drop = FALSE]
mat_thin800K <- mat_thin800K[rownames(mat_thin800K) %in% LH_Data$ID, , drop = FALSE]
mat_thin900K <- mat_thin900K[rownames(mat_thin900K) %in% LH_Data$ID, , drop = FALSE]

# 208 indivs
dim(mat_thin25K) # 1100
dim(mat_thin50K) # 1041
dim(mat_thin75K) # 1005
dim(mat_thin100K) # 943
dim(mat_thin200K) # 808
dim(mat_thin300K) # 711
dim(mat_thin400K) # 648
dim(mat_thin500K) # 594
dim(mat_thin600K) # 546
dim(mat_thin700K) # 509
dim(mat_thin800K) # 476
dim(mat_thin900K) # 446





##### ---- Final GenoMat Maintenence and Checks ---- #####

# check for all heterozygous sites since those likely will not be informative.
all_ones_25 <- apply(mat_thin25K, 2, function(col) all(col == 1))
all_ones_50 <- apply(mat_thin50K, 2, function(col) all(col == 1))
all_ones_75 <- apply(mat_thin75K, 2, function(col) all(col == 1))
all_ones_1 <- apply(mat_thin100K, 2, function(col) all(col == 1))
all_ones_2 <- apply(mat_thin200K, 2, function(col) all(col == 1))
all_ones_3 <- apply(mat_thin300K, 2, function(col) all(col == 1))
all_ones_4 <- apply(mat_thin400K, 2, function(col) all(col == 1))
all_ones_5 <- apply(mat_thin500K, 2, function(col) all(col == 1))
all_ones_6 <- apply(mat_thin600K, 2, function(col) all(col == 1))
all_ones_7 <- apply(mat_thin700K, 2, function(col) all(col == 1))
all_ones_8 <- apply(mat_thin800K, 2, function(col) all(col == 1))
all_ones_9 <- apply(mat_thin900K, 2, function(col) all(col == 1))

# remove all het. sites
mat_thin25K <- mat_thin25K[, !all_ones_25]
mat_thin50K <- mat_thin50K[, !all_ones_50]
mat_thin75K <- mat_thin75K[, !all_ones_75]
mat_thin100K <- mat_thin100K[, !all_ones_1]
mat_thin200K <- mat_thin200K[, !all_ones_2]
mat_thin300K <- mat_thin300K[, !all_ones_3]
mat_thin400K <- mat_thin400K[, !all_ones_4]
mat_thin500K <- mat_thin500K[, !all_ones_5]
mat_thin600K <- mat_thin600K[, !all_ones_6]
mat_thin700K <- mat_thin700K[, !all_ones_7]
mat_thin800K <- mat_thin800K[, !all_ones_8]
mat_thin900K <- mat_thin900K[, !all_ones_9]

dim(mat_thin25K) # 1100 <- 1099
dim(mat_thin50K) # 1041 <- 1040
dim(mat_thin75K) # 1005 <- 1004
dim(mat_thin100K) # 943 <- 942
dim(mat_thin200K) # 808 <- 807
dim(mat_thin300K) # 711 <- 710
dim(mat_thin400K) # 648 <- 647
dim(mat_thin500K) # 594 <- 593
dim(mat_thin600K) # 546 <- 545
dim(mat_thin700K) # 509 <- 508
dim(mat_thin800K) # 476 <- 476
dim(mat_thin900K) # 446 <- 445
# lost one site in each, except for 800K.

check_thin25K <- CheckGeno(mat_thin25K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin50K <- CheckGeno(mat_thin50K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin75K <- CheckGeno(mat_thin75K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin100K <- CheckGeno(mat_thin100K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin200K <- CheckGeno(mat_thin200K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin300K <- CheckGeno(mat_thin300K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin400K <- CheckGeno(mat_thin400K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin500K <- CheckGeno(mat_thin500K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin600K <- CheckGeno(mat_thin600K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin700K <- CheckGeno(mat_thin700K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin800K <- CheckGeno(mat_thin800K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin900K <- CheckGeno(mat_thin900K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))





















##### ---- Running GetMaybeRel(), Checking with cross_lookup ---- #####

# info here: alright, when i did this for HOURS on Sequoia_ContamFilt_mindep8.R,
# i had so many runs of this, and that script got up to over 2000 lines. since i
# am already planning to store all of this information in sequoia_params.xlsx, i
# am going to just use this one template and sub out the 1-9 matrices here. see
# results in that table.

# set error rate
erm1 <- ErrToM(Err = 0.01)
erm25 <- ErrToM(Err = 0.025)
erm5 <- ErrToM(Err = 0.05)

# run gmr
gmr_thin900K <- GetMaybeRel(GenoM = check_thin900K,                             # here (2)
                               Err = erm5,                                      # error
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_thin900K))               # here (1)

#unique test f1's in the trios:
length(unique(gmr_thin900K[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_thin900K[["MaybeTrio"]]$id)])) # here (2)

# how many unique focal indivs are in the trios?
length(unique(gmr_thin900K[["MaybeTrio"]]$id))                                  # here (1)
# see who it is
table(unique(gmr_thin900K[["MaybeTrio"]]$id))                                   # here (1)

trios_thin900K <- gmr_thin900K[["MaybeTrio"]]                                   # here (2)
head(trios_thin900K)                                                            # here (1)

# read in lookup table
cross_lookup <- read.csv(file = "true_cross_lookup.csv", header = TRUE)
cross_lookup <- cross_lookup[,-1]

# this is a ridiculously effective set of piped functions here:
trios_thin900K <- trios_thin900K %>%                                            # here (2)
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

head(trios_thin900K)                                                            # here (1)
dim(trios_thin900K)                                                             # here (1)
table(trios_thin900K$valid_cross)                                               # here (1)

trios_thin900K_checked <- trios_thin900K %>%                                    # here (2)
  select(id, parent1, parent2, pair, valid_cross, everything())

trios_thin900K_checked <- trios_thin900K_checked %>%                            # here (2)
  filter(!LLRparent1 %in% c(555, -555),
         !LLRparent2 %in% c(555, -555),
         !LLRpair    %in% c(555, -555))
dim(trios_thin900K_checked)                                                     # here (1)
table(trios_thin900K_checked$valid_cross)                                       # here (1)

length(unique(trios_thin900K_checked$id[grepl("^SAR_15_67", trios_thin900K_checked$id)])) # here (2)












