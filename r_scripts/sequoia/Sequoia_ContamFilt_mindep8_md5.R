## Sequoia_ContamFilt_mindep8_md5.R by SPJ 091225
## PURPOSE: I have compiled data of sequoia results for various n snps, md, and
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

# NOTES 092525: one concern that's arisen are loci that are out of HWE. seems easy
# enough for me to calculate the HWE.p values from the HWE chisq test and see 1.
# how many loci are OUT of HWE, and 2. the distrubution of scores for each datset.
# for this, i need all 17 datasets that have md 5%. if you see in the notes above,
# the thinned 1M to 5M were done in another script, Sequoia_ContamFilt_mindep8.R,
# because those were the first ones I tested for the sequoia_params table, then I
# did all the md10 ones to get more SNPs, then came back and did all the rest of the
# md5 ones in here. I'm going to go to Sequoia_ContamFilt_mindep8.R, grab the code to
# construct the matrices and get them ready for GetMaybeRel(), then I'm going to run
# SnpStats on them to get the information that I'm looking for.

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

# these are read in from Sequoia_ContamFilt_mindep8.R, names changed from including min8 to not.

mat_thin1M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin1M.012_conv",
                              header = FALSE, sep = "\t", na.strings = c("NA", "-1")) 

mat_thin2M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin2M.012_conv",
                              header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin3M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin3M.012_conv",
                              header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin4M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin4M.012_conv",
                              header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin5M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin5M.012_conv",
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
mat_thin1M <- mat_thin1M[, -1]
mat_thin2M <- mat_thin2M[, -1]
mat_thin3M <- mat_thin3M[, -1]
mat_thin4M <- mat_thin4M[, -1]
mat_thin5M <- mat_thin5M[, -1]

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
mat_thin1M <- as.matrix(mat_thin1M)
mat_thin2M <- as.matrix(mat_thin2M)
mat_thin3M <- as.matrix(mat_thin3M)
mat_thin4M <- as.matrix(mat_thin4M)
mat_thin5M <- as.matrix(mat_thin5M)

dim(mat_thin25K)
dim(mat_thin50K)
dim(mat_thin75K)
dim(mat_thin100K)
dim(mat_thin200K)
dim(mat_thin300K)
dim(mat_thin400K)
dim(mat_thin500K)
dim(mat_thin600K)
dim(mat_thin700K)
dim(mat_thin800K)
dim(mat_thin900K)
dim(mat_thin1M)
dim(mat_thin2M)
dim(mat_thin3M)
dim(mat_thin4M)
dim(mat_thin5M)


# yes, i am aware that these are all the same. trying to be thorough. would hate 
# to mess up a matrix because i was lazy on a step like this.
ind25K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin25K.012.indv", header = FALSE)
ind50K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin50K.012.indv", header = FALSE)
ind75K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin75K.012.indv", header = FALSE)
ind100K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012.indv", header = FALSE)
ind200K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin200K.012.indv", header = FALSE)
ind300K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin300K.012.indv", header = FALSE)
ind400K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin400K.012.indv", header = FALSE)
ind500K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin500K.012.indv", header = FALSE)
ind600K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin600K.012.indv", header = FALSE)
ind700K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin700K.012.indv", header = FALSE)
ind800K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin800K.012.indv", header = FALSE)
ind900K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin900K.012.indv", header = FALSE)
ind1M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin1M.012.indv", header = FALSE)
ind2M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin2M.012.indv", header = FALSE)
ind3M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin3M.012.indv", header = FALSE)
ind4M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin4M.012.indv", header = FALSE)
ind5M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin5M.012.indv", header = FALSE)

# rename first column to something meaningful
ind25K <- ind25K %>% 
  rename(sample = V1)
ind50K <- ind50K %>% 
  rename(sample = V1)
ind75K <- ind75K %>% 
  rename(sample = V1)
ind100K <-  ind100K %>% 
  rename(sample = V1)
ind200K <-  ind200K %>% 
  rename(sample = V1)
ind300K <-  ind300K %>% 
  rename(sample = V1)
ind400K <-  ind400K %>% 
  rename(sample = V1)
ind500K <-  ind500K %>% 
  rename(sample = V1)
ind600K <-  ind600K %>% 
  rename(sample = V1)
ind700K <-  ind700K %>% 
  rename(sample = V1)
ind800K <-  ind800K %>% 
  rename(sample = V1)
ind900K <-  ind900K %>% 
  rename(sample = V1)
ind1M <- ind1M %>% 
  rename(sample = V1)
ind2M <- ind2M %>% 
  rename(sample = V1)
ind3M <- ind3M %>% 
  rename(sample = V1)
ind4M <- ind4M %>% 
  rename(sample = V1)
ind5M <- ind5M %>% 
  rename(sample = V1)

# establish sample identities in the geno_mat
rownames(mat_thin25K) <- ind25K$sample
rownames(mat_thin50K) <- ind50K$sample
rownames(mat_thin75K) <- ind75K$sample
rownames(mat_thin100K) <- ind100K$sample
rownames(mat_thin200K) <- ind200K$sample
rownames(mat_thin300K) <- ind300K$sample
rownames(mat_thin400K) <- ind400K$sample
rownames(mat_thin500K) <- ind500K$sample
rownames(mat_thin600K) <- ind600K$sample
rownames(mat_thin700K) <- ind700K$sample
rownames(mat_thin800K) <- ind800K$sample
rownames(mat_thin900K) <- ind900K$sample
rownames(mat_thin1M) <- ind1M$sample
rownames(mat_thin2M) <- ind2M$sample
rownames(mat_thin3M) <- ind3M$sample
rownames(mat_thin4M) <- ind4M$sample
rownames(mat_thin5M) <- ind5M$sample

# read in scaffolds and positions
pos25K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin25K.012.pos", header = FALSE)
pos50K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin50K.012.pos", header = FALSE)
pos75K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin75K.012.pos", header = FALSE)
pos100K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012.pos", header = FALSE)
pos200K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin200K.012.pos", header = FALSE)
pos300K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin300K.012.pos", header = FALSE)
pos400K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin400K.012.pos", header = FALSE)
pos500K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin500K.012.pos", header = FALSE)
pos600K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin600K.012.pos", header = FALSE)
pos700K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin700K.012.pos", header = FALSE)
pos800K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin800K.012.pos", header = FALSE)
pos900K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin900K.012.pos", header = FALSE)
pos1M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin1M.012.pos", header = FALSE)
pos2M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin2M.012.pos", header = FALSE)
pos3M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin3M.012.pos", header = FALSE)
pos4M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin4M.012.pos", header = FALSE)
pos5M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin5M.012.pos", header = FALSE)

# create full positions by combining the two columns
pos25K$position <- paste(pos25K$V1, pos25K$V2, sep = "_")
pos50K$position <- paste(pos50K$V1, pos50K$V2, sep = "_")
pos75K$position <- paste(pos75K$V1, pos75K$V2, sep = "_")
pos100K$position <-  paste(pos100K$V1, pos100K$V2, sep = "_")
pos200K$position <-  paste(pos200K$V1, pos200K$V2, sep = "_")
pos300K$position <-  paste(pos300K$V1, pos300K$V2, sep = "_")
pos400K$position <-  paste(pos400K$V1, pos400K$V2, sep = "_")
pos500K$position <-  paste(pos500K$V1, pos500K$V2, sep = "_")
pos600K$position <-  paste(pos600K$V1, pos600K$V2, sep = "_")
pos700K$position <-  paste(pos700K$V1, pos700K$V2, sep = "_")
pos800K$position <-  paste(pos800K$V1, pos800K$V2, sep = "_")
pos900K$position <-  paste(pos900K$V1, pos900K$V2, sep = "_")
pos1M$position <- paste(pos1M$V1, pos1M$V2, sep = "_")
pos2M$position <- paste(pos2M$V1, pos2M$V2, sep = "_")
pos3M$position <- paste(pos3M$V1, pos3M$V2, sep = "_")
pos4M$position <- paste(pos4M$V1, pos4M$V2, sep = "_")
pos5M$position <- paste(pos5M$V1, pos5M$V2, sep = "_")

# set positions as column names of geno_mat
colnames(mat_thin25K) <- pos25K$position
colnames(mat_thin50K) <- pos50K$position
colnames(mat_thin75K) <- pos75K$position
colnames(mat_thin100K) <- pos100K$position
colnames(mat_thin200K) <- pos200K$position
colnames(mat_thin300K) <- pos300K$position
colnames(mat_thin400K) <- pos400K$position
colnames(mat_thin500K) <- pos500K$position
colnames(mat_thin600K) <- pos600K$position
colnames(mat_thin700K) <- pos700K$position
colnames(mat_thin800K) <- pos800K$position
colnames(mat_thin900K) <- pos900K$position
colnames(mat_thin1M) <- pos1M$position
colnames(mat_thin2M) <- pos2M$position
colnames(mat_thin3M) <- pos3M$position
colnames(mat_thin4M) <- pos4M$position
colnames(mat_thin5M) <- pos5M$position





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
LH_Data$BirthYear[1:95] <- 1
LH_Data$BirthYear[96:nrow(LH_Data)] <- 0
# THIS IS READY NOW FOR GETMAYBEREL

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
mat_thin1M <- mat_thin1M[rownames(mat_thin1M) %in% LH_Data$ID, , drop = FALSE]
mat_thin2M <- mat_thin2M[rownames(mat_thin2M) %in% LH_Data$ID, , drop = FALSE]
mat_thin3M <- mat_thin3M[rownames(mat_thin3M) %in% LH_Data$ID, , drop = FALSE]
mat_thin4M <- mat_thin4M[rownames(mat_thin4M) %in% LH_Data$ID, , drop = FALSE]
mat_thin5M <- mat_thin5M[rownames(mat_thin5M) %in% LH_Data$ID, , drop = FALSE]

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
dim(mat_thin1M) # 418
dim(mat_thin2M) # 273
dim(mat_thin3M) # 211
dim(mat_thin4M) # 175
dim(mat_thin5M) # 151



##### ---- HWE Checks and Final GenoMat Maintenence ---- #####

# check for all heterozygous sites since those likely will not be informative.
all_ones_25 <- apply(mat_thin25K, 2, function(col) all(col == 1))
all_ones_50 <- apply(mat_thin50K, 2, function(col) all(col == 1))
all_ones_75 <- apply(mat_thin75K, 2, function(col) all(col == 1))
all_ones_100 <- apply(mat_thin100K, 2, function(col) all(col == 1))
all_ones_200 <- apply(mat_thin200K, 2, function(col) all(col == 1))
all_ones_300 <- apply(mat_thin300K, 2, function(col) all(col == 1))
all_ones_400 <- apply(mat_thin400K, 2, function(col) all(col == 1))
all_ones_500 <- apply(mat_thin500K, 2, function(col) all(col == 1))
all_ones_600 <- apply(mat_thin600K, 2, function(col) all(col == 1))
all_ones_700 <- apply(mat_thin700K, 2, function(col) all(col == 1))
all_ones_800 <- apply(mat_thin800K, 2, function(col) all(col == 1))
all_ones_900 <- apply(mat_thin900K, 2, function(col) all(col == 1))
all_ones_1 <- apply(mat_thin1M, 2, function(col) all(col == 1))
all_ones_2 <- apply(mat_thin2M, 2, function(col) all(col == 1))
all_ones_3 <- apply(mat_thin3M, 2, function(col) all(col == 1))
all_ones_4 <- apply(mat_thin4M, 2, function(col) all(col == 1))
all_ones_5 <- apply(mat_thin5M, 2, function(col) all(col == 1))

# remove all het. sites
mat_thin25K <- mat_thin25K[, !all_ones_25]
mat_thin50K <- mat_thin50K[, !all_ones_50]
mat_thin75K <- mat_thin75K[, !all_ones_75]
mat_thin100K <- mat_thin100K[, !all_ones_100]
mat_thin200K <- mat_thin200K[, !all_ones_200]
mat_thin300K <- mat_thin300K[, !all_ones_300]
mat_thin400K <- mat_thin400K[, !all_ones_400]
mat_thin500K <- mat_thin500K[, !all_ones_500]
mat_thin600K <- mat_thin600K[, !all_ones_600]
mat_thin700K <- mat_thin700K[, !all_ones_700]
mat_thin800K <- mat_thin800K[, !all_ones_800]
mat_thin900K <- mat_thin900K[, !all_ones_900]
mat_thin1M <- mat_thin1M[, !all_ones_1]
mat_thin2M <- mat_thin2M[, !all_ones_2]
mat_thin3M <- mat_thin3M[, !all_ones_3] 
mat_thin4M <- mat_thin4M[, !all_ones_4]
mat_thin5M <- mat_thin5M[, !all_ones_5] 
  
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
dim(mat_thin1M) # 418 <- 418
dim(mat_thin2M) # 273 <- 272
dim(mat_thin3M) # 211 <- 211
dim(mat_thin4M) # 175 <- 175
dim(mat_thin5M) # 151 <- 151
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
check_thin1M <- CheckGeno(mat_thin1M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin2M <- CheckGeno(mat_thin2M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin3M <- CheckGeno(mat_thin3M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin4M <- CheckGeno(mat_thin4M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin5M <- CheckGeno(mat_thin5M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))


##### ---- SNP Stats, checking for HWE, and filtering for HWE
stats_25K <- data.frame(SnpStats(GenoM = mat_thin25K, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_50K <- data.frame(SnpStats(GenoM = mat_thin50K, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_75K <- data.frame(SnpStats(GenoM = mat_thin75K, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_100K <- data.frame(SnpStats(GenoM = mat_thin100K, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_200K <- data.frame(SnpStats(GenoM = mat_thin200K, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_300K <- data.frame(SnpStats(GenoM = mat_thin300K, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_400K <- data.frame(SnpStats(GenoM = mat_thin400K, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_500K <- data.frame(SnpStats(GenoM = mat_thin500K, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_600K <- data.frame(SnpStats(GenoM = mat_thin600K, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_700K <- data.frame(SnpStats(GenoM = mat_thin700K, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_800K <- data.frame(SnpStats(GenoM = mat_thin800K, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_900K <- data.frame(SnpStats(GenoM = mat_thin900K, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_1M <- data.frame(SnpStats(GenoM = mat_thin1M, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_2M <- data.frame(SnpStats(GenoM = mat_thin2M, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_3M <- data.frame(SnpStats(GenoM = mat_thin3M, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_4M <- data.frame(SnpStats(GenoM = mat_thin4M, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))
stats_5M <- data.frame(SnpStats(GenoM = mat_thin5M, Pedigree = NULL, Plot = TRUE, quiet = FALSE, calc_HWE = TRUE))

# create a list of each of the objects
GenoM_list <- list(
  stats_25K = stats_25K, 
  stats_50K = stats_50K, 
  stats_75K = stats_75K, 
  stats_100K = stats_100K, 
  stats_200K = stats_200K, 
  stats_300K = stats_300K, 
  stats_400K = stats_400K, 
  stats_500K = stats_500K, 
  stats_600K = stats_600K, 
  stats_700K = stats_700K, 
  stats_800K = stats_800K, 
  stats_900K = stats_900K,
  stats_1M = stats_1M,
  stats_2M = stats_2M,
  stats_3M = stats_3M,
  stats_4M = stats_4M,
  stats_5M = stats_5M
)

# create a summary list, that applies over each GenoM, gives desired info for each
summary_list <- lapply(GenoM_list, function(df) {
  n_snps <- nrow(df)
  n_out <- sum(df$HWE.p < 0.05, na.rm = TRUE)
  prop_out <- n_out / n_snps
  data.frame(
    n_snps = n_snps,
    n_out_HWE = n_out,
    prop_out_HWE = prop_out
  )
})

# then rbind them together
summary_df <- do.call(rbind, summary_list)

# add an identifier column and place it in front
summary_df$GenoM <- names(GenoM_list)
summary_df <- summary_df %>% 
    select(GenoM, everything())
summary_df # check

# plot distributions of HWE p values for each GenoM.
library(dplyr)
library(ggplot2)

hwe_df <- bind_rows(GenoM_list, .id = "GenoM") %>% 
  mutate(minus_log10_p = log10(HWE.p))

ggplot(hwe_df, aes(x = GenoM, y = minus_log10_p)) +
  geom_violin(fill = "skyblue", alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed", color = "red") +
  labs(
    title = "Dist. of HWE p-values across datasets",
    x = "Dataset",
    y = expression(log[10](HWE~p))
  ) +
  theme_bw()

# alright, now we want to filter out the SNPs in each dataset that are out of HWE.
hwe_threshold <- 0.05 # set threshold

# create a list of the unique character strings for the list
datasets <- c("25K", "50K", "75K", "100K", "200K", "300K", "400K", "500K", "600K", 
              "700K", "800K", "900K", "1M", "2M", "3M", "4M", "5M")

# for each of these unique character strings (thinning distances)
for(d in datasets){
  stats_obj <- get(paste0("stats_", d)) # get the snpstats object that corresponds to that distance
  mat_obj <- get(paste0("mat_thin", d)) # and get the GenoM that corresponds to that distance
  
  # create a vector of the rownames of the stats object, where the HWE.p value is below 0.05
  snps_fail_hwe <- rownames(stats_obj)[stats_obj$HWE.p < hwe_threshold]
  
  # creat a filtered matirx, from the original matrix object, keep all rows, but remove columns
  # where the column name is a character string that's in the snps_fail_hwe object
  mat_filtered <- mat_obj[, !colnames(mat_obj) %in% snps_fail_hwe]
  
  # now assign the filtered matrix to a new object that's named appropriately, with a _HWEfilt extension
  assign(paste0("mat_thin", d, "_HWEfilt"), mat_filtered)
}

# now let's recreate the check_thin objects FROM the HWE filtered objects
check_thin25K <- CheckGeno(mat_thin25K_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin50K <- CheckGeno(mat_thin50K_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin75K <- CheckGeno(mat_thin75K_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin100K <- CheckGeno(mat_thin100K_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin200K <- CheckGeno(mat_thin200K_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin300K <- CheckGeno(mat_thin300K_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin400K <- CheckGeno(mat_thin400K_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin500K <- CheckGeno(mat_thin500K_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin600K <- CheckGeno(mat_thin600K_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin700K <- CheckGeno(mat_thin700K_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin800K <- CheckGeno(mat_thin800K_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin900K <- CheckGeno(mat_thin900K_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin1M <- CheckGeno(mat_thin1M_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin2M <- CheckGeno(mat_thin2M_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin3M <- CheckGeno(mat_thin3M_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin4M <- CheckGeno(mat_thin4M_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_thin5M <- CheckGeno(mat_thin5M_HWEfilt, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))

dim(check_thin25K) # 208 733
dim(check_thin50K) # 208 694
dim(check_thin75K) # 208 676
dim(check_thin100K) # 208 633
dim(check_thin200K) # 208 541
dim(check_thin300K) # 208 482
dim(check_thin400K) # 208 445
dim(check_thin500K) # 208 406
dim(check_thin600K) # 208 368
dim(check_thin700K) # 208 347
dim(check_thin800K) # 208 331
dim(check_thin900K) # 208 300
dim(check_thin1M) # 208 288
dim(check_thin2M) # 208 183
dim(check_thin3M) # 208 133
dim(check_thin4M) # 208 106
dim(check_thin5M) # 208 101


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
gmr_thin4M <- GetMaybeRel(GenoM = check_thin4M,                             # here (2)
                               Err = erm25,                                      # error
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               # Herm = "no",
                               MaxPairs = 7*nrow(check_thin4M))               # here (1)

#unique test f1's in the trios:
length(unique(gmr_thin4M[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_thin4M[["MaybeTrio"]]$id)])) # here (2)

# how many unique focal indivs are in the trios?
length(unique(gmr_thin4M[["MaybeTrio"]]$id))                                  # here (1)
# see who it is
table(unique(gmr_thin4M[["MaybeTrio"]]$id))                                   # here (1)

trios_thin4M <- gmr_thin4M[["MaybeTrio"]]                                   # here (2)
head(trios_thin4M)                                                            # here (1)

# read in lookup table
cross_lookup <- read.csv(file = "true_cross_lookup.csv", header = TRUE)
cross_lookup <- cross_lookup[,-1]

# this is a ridiculously effective set of piped functions here:
trios_thin4M <- trios_thin4M %>%                                            # here (2)
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

head(trios_thin4M)                                                            # here (1)
dim(trios_thin4M)                                                             # here (1)
table(trios_thin4M$valid_cross)                                               # here (1)

trios_thin4M_checked <- trios_thin4M %>%                                    # here (2)
  select(id, parent1, parent2, pair, valid_cross, everything())

trios_thin4M_checked <- trios_thin4M_checked %>%                            # here (2)
  filter(!LLRparent1 %in% c(555, -555),
         !LLRparent2 %in% c(555, -555),
         !LLRpair    %in% c(555, -555))
dim(trios_thin4M_checked)                                                     # here (1)
table(trios_thin4M_checked$valid_cross)                                       # here (1)

length(unique(trios_thin4M_checked$id[grepl("^SAR_15_67", trios_thin4M_checked$id)])) # here (2)








