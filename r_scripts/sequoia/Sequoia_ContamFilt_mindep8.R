## Sequoia_ContamFilt_mindep8.R by SPJ 081325
## PURPOSE: to implement suggestions from Jisca into a new dataset that has been
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
# q40 - keep sites with quality > 40 (mean base quality? or is this mean q for the read?)
# mindep8 - minimum mean read depth 4 or above (removes low quality sites, low confidence calls)
# maxdep75 - maximum mean read depth 75 or lower (removes paralogs)
# maf30 - minor allele frequency 30% or higher
# miss95 - each snp needs to have 5% or less missing data

# and thinned to various degrees to filter for LD

# these vcf -> .012 came from me wanting to be more confident in my genotype calls
# which is why q went up from 20 to 40, and mindep went up to 8.

# before 08/22/25 is me prepping the matrices and running them on sequoia's 
# latest version on CRAN at the time of 08/2025 (version 2.11.2)

# starting on 08/22/25 is where I start experimenting with 
# older versions (2.8.3) to see if that will help increase the number of assignments
# that I'm able to make.

# note: i had to restart R studio, prior to that point, so I first had to rebuild
# the matrices in sequoia 2.8.3 but that'll only affect sequoia functions such as
# CheckGeno() and GetMaybeRel().

# then on 08/25/25 i started working with the newest version (3.0.3)

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

mat_min8_thin1M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin1M.012_conv",
                              header = FALSE, sep = "\t", na.strings = c("NA", "-1")) 

mat_min8_thin2M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin2M.012_conv",
                              header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_min8_thin3M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin3M.012_conv",
                              header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_min8_thin4M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin4M.012_conv",
                              header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_min8_thin5M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin5M.012_conv",
                              header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

# remove first column. ascending numbers 1 - 1184
mat_min8_thin1M <- mat_min8_thin1M[, -1]
mat_min8_thin2M <- mat_min8_thin2M[, -1]
mat_min8_thin3M <- mat_min8_thin3M[, -1]
mat_min8_thin4M <- mat_min8_thin4M[, -1]
mat_min8_thin5M <- mat_min8_thin5M[, -1]

# turn them to matrices
mat_min8_thin1M <- as.matrix(mat_min8_thin1M)
mat_min8_thin2M <- as.matrix(mat_min8_thin2M)
mat_min8_thin3M <- as.matrix(mat_min8_thin3M)
mat_min8_thin4M <- as.matrix(mat_min8_thin4M)
mat_min8_thin5M <- as.matrix(mat_min8_thin5M)

# yes, i am aware that these are all the same. trying to be thorough. would hate 
# to mess up a matrix because i was lazy on a step like this.
ind_min8_thin1M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin1M.012.indv", header = FALSE)
ind_min8_thin2M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin2M.012.indv", header = FALSE)
ind_min8_thin3M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin3M.012.indv", header = FALSE)
ind_min8_thin4M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin4M.012.indv", header = FALSE)
ind_min8_thin5M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin5M.012.indv", header = FALSE)

# rename first column to something meaningful
ind_min8_thin1M <- ind_min8_thin1M  %>% 
  rename(sample = V1)
ind_min8_thin2M <- ind_min8_thin2M  %>% 
  rename(sample = V1)
ind_min8_thin3M <- ind_min8_thin3M  %>% 
  rename(sample = V1)
ind_min8_thin4M <- ind_min8_thin4M  %>% 
  rename(sample = V1)
ind_min8_thin5M <- ind_min8_thin5M  %>% 
  rename(sample = V1)

# establish sample identities in the geno_mat
rownames(mat_min8_thin1M) <- ind_min8_thin1M$sample
rownames(mat_min8_thin2M) <- ind_min8_thin2M$sample
rownames(mat_min8_thin3M) <- ind_min8_thin3M$sample
rownames(mat_min8_thin4M) <- ind_min8_thin4M$sample
rownames(mat_min8_thin5M) <- ind_min8_thin5M$sample

# read in scaffolds and positions
pos_min8_thin1M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin1M.012.pos", header = FALSE)
pos_min8_thin2M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin2M.012.pos", header = FALSE)
pos_min8_thin3M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin3M.012.pos", header = FALSE)
pos_min8_thin4M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin4M.012.pos", header = FALSE)
pos_min8_thin5M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin5M.012.pos", header = FALSE)

# create full positions by combining the two columns
pos_min8_thin1M$position <- paste(pos_min8_thin1M$V1, pos_min8_thin1M$V2, sep = "_")
pos_min8_thin2M$position <- paste(pos_min8_thin2M$V1, pos_min8_thin2M$V2, sep = "_")
pos_min8_thin3M$position <- paste(pos_min8_thin3M$V1, pos_min8_thin3M$V2, sep = "_")
pos_min8_thin4M$position <- paste(pos_min8_thin4M$V1, pos_min8_thin4M$V2, sep = "_")
pos_min8_thin5M$position <- paste(pos_min8_thin5M$V1, pos_min8_thin5M$V2, sep = "_")

# set positions as column names of geno_mat
colnames(mat_min8_thin1M) <- pos_min8_thin1M$position
colnames(mat_min8_thin2M) <- pos_min8_thin2M$position
colnames(mat_min8_thin3M) <- pos_min8_thin3M$position
colnames(mat_min8_thin4M) <- pos_min8_thin4M$position
colnames(mat_min8_thin5M) <- pos_min8_thin5M$position

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
mat_min8_thin1M <- mat_min8_thin1M[rownames(mat_min8_thin1M) %in% LH_Data$ID, , drop = FALSE]
dim(mat_min8_thin1M) # 208 indivs, 418 loci

mat_min8_thin2M <- mat_min8_thin2M[rownames(mat_min8_thin2M) %in% LH_Data$ID, , drop = FALSE]
dim(mat_min8_thin2M) # 208 indivs, 273 loci

mat_min8_thin3M <- mat_min8_thin3M[rownames(mat_min8_thin3M) %in% LH_Data$ID, , drop = FALSE]
dim(mat_min8_thin3M) # 208 indivs, 211 loci

mat_min8_thin4M <- mat_min8_thin4M[rownames(mat_min8_thin4M) %in% LH_Data$ID, , drop = FALSE]
dim(mat_min8_thin4M) # 208 indivs, 175 loci

mat_min8_thin5M <- mat_min8_thin5M[rownames(mat_min8_thin5M) %in% LH_Data$ID, , drop = FALSE]
dim(mat_min8_thin5M) # 208 indivs, 151 loci

# check for all heterozygous sites since those likely will not be informative.
all_ones_min8_thin1M <- apply(mat_min8_thin1M, 2, function(col) all(col == 1))
all_ones_min8_thin2M <- apply(mat_min8_thin2M, 2, function(col) all(col == 1))
all_ones_min8_thin3M <- apply(mat_min8_thin3M, 2, function(col) all(col == 1))
all_ones_min8_thin4M <- apply(mat_min8_thin4M, 2, function(col) all(col == 1))
all_ones_min8_thin5M <- apply(mat_min8_thin5M, 2, function(col) all(col == 1))

mat_min8_thin1M <- mat_min8_thin1M[, !all_ones_min8_thin1M]
dim(mat_min8_thin1M) # none lost, 418

mat_min8_thin2M <- mat_min8_thin2M[, !all_ones_min8_thin2M]
dim(mat_min8_thin2M) # ONE lost, 272

mat_min8_thin3M <- mat_min8_thin3M[, !all_ones_min8_thin3M]
dim(mat_min8_thin3M) # none lost, 211

mat_min8_thin4M <- mat_min8_thin4M[, !all_ones_min8_thin4M]
dim(mat_min8_thin4M) # none lost, 175

mat_min8_thin5M <- mat_min8_thin5M[, !all_ones_min8_thin5M]
dim(mat_min8_thin5M) # none lost, 151

# check genotype matrix for samples/loci to be excluded
check_min8_thin1M <- CheckGeno(mat_min8_thin1M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_min8_thin2M <- CheckGeno(mat_min8_thin2M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_min8_thin3M <- CheckGeno(mat_min8_thin3M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_min8_thin4M <- CheckGeno(mat_min8_thin4M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check_min8_thin5M <- CheckGeno(mat_min8_thin5M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))

# alright, let's now remove any sites with missing data and see how many are retained
# create an index of sites where there are no -9 values
completesites_min8_thin1M <- apply(check_min8_thin1M, 2, function(col) all(col != -9))
completesites_min8_thin2M <- apply(check_min8_thin2M, 2, function(col) all(col != -9))
completesites_min8_thin3M <- apply(check_min8_thin3M, 2, function(col) all(col != -9))
completesites_min8_thin4M <- apply(check_min8_thin4M, 2, function(col) all(col != -9))
completesites_min8_thin5M <- apply(check_min8_thin5M, 2, function(col) all(col != -9))

# subset the matrix to keep only those columns
check_min8_thin1M_nomiss <- check_min8_thin1M[, completesites_min8_thin1M]
dim(check_min8_thin1M_nomiss) # 137 SNPs

check_min8_thin2M_nomiss <- check_min8_thin2M[, completesites_min8_thin2M]
dim(check_min8_thin2M_nomiss) # 87 SNPs

check_min8_thin3M_nomiss <- check_min8_thin3M[, completesites_min8_thin3M]
dim(check_min8_thin3M_nomiss) # 79 SNPs

check_min8_thin4M_nomiss <- check_min8_thin4M[, completesites_min8_thin4M]
dim(check_min8_thin4M_nomiss) # 53 SNPs

check_min8_thin5M_nomiss <- check_min8_thin5M[, completesites_min8_thin5M]
dim(check_min8_thin5M_nomiss) # 47 SNPs

# let's first try all of these with errors of 1, 2.5, and 5%
# ten matrices with these three error combinations = 30 runs.
error_rate <- c(0.015, 0.015, 0.015)

gmr_min8_thin1M_nomiss <- GetMaybeRel(GenoM = check_min8_thin1M_nomiss,
                            Err = error_rate,
                            # SeqList = outfull,
                            Module = "par",
                            # MaxMismatch = NA,
                            Complex = "simp",
                            LifeHistData = LH_Data,
                            quiet = FALSE,
                            Tassign = 1.0,
                            # Tfilter = -100,
                            MaxPairs = 7*nrow(check_min8_thin1M_nomiss))
# Error 1%: 165 PO pairs, 32 others, 83 trios. 137 SNPs.
# Error 1.5%: 165 PO pairs, 37 others, 82 trios. 137 SNPs.
# Error 2.5%: 175 PO pairs, 52 others, 90 trios. 137 SNPs.
# Error 3%: 155 PO pairs, 61 others, 93 trios. 137 SNPs.
# Error 5%: 163 PO pairs, 91 others, 113 trios. 137 SNPs.

gmr_min8_thin2M_nomiss <- GetMaybeRel(GenoM = check_min8_thin2M_nomiss,
                                      Err = error_rate,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin2M_nomiss))
# Error 1%: 158 PO pairs, 57 others, 81 trios. 87 SNPs.
# Error 1.5%: 168 PO pairs, 65 others, 89 trios. 87 SNPs.
# Error 2.5%: 177 PO pairs, 80 others, 101 trios. 87 SNPs.
# Error 3%: 174 PO pairs, 95 other,s 100 trios. 87 SNPs.
# Error 5%: 168 PO pairs, 135 others, 127 trios. 87 SNPs.

gmr_min8_thin3M_nomiss <- GetMaybeRel(GenoM = check_min8_thin3M_nomiss,
                                      Err = error_rate,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin3M_nomiss))
# Error 1%: 149 PO pairs, 65 others, 78 trios. 79 SNPs.
# Error 1.5%: 158 PO pairs, 75 others, 85 trios. 79 SNPs.
# Error 2.5%: 162 PO pairs, 95 others, 101 trios. 79 SNPs.
# Error 3%: 161 PO pairs, 105 other,s 102 trios. 79 SNPs.
# Error 5%: 155 PO pairs, 138 others, 133 trios. 79 SNPs.

gmr_min8_thin4M_nomiss <- GetMaybeRel(GenoM = check_min8_thin4M_nomiss,
                                      Err = error_rate,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin4M_nomiss))
# Error 1%: 186 PO pairs, 81 others, 84 trios. 53 SNPs.
# Error 1.5%: 190 PO pairs, 87 others, 99 trios. 53 SNPs.
# Error 2.5%: 189 PO pairs, 109 others, 119 trios. 53 SNPs.
# Error 3%: 183 PO pairs, 118 others, 127 trios. 53 SNPs.
# Error 5%: 160 PO pairs, 129 others, 134 trios. 53 SNPs.

gmr_min8_thin5M_nomiss <- GetMaybeRel(GenoM = check_min8_thin5M_nomiss,
                                      Err = error_rate,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin5M_nomiss))
# Error 1%: 186 PO pairs, 81 others, 72 trios. 47 SNPs.
# Error 1.5%: 181 PO pairs, 92 others, 74 trios. 47 SNPs.
# Error 2.5%: 168 PO pairs, 102 others, 78 trios. 47 SNPs.
# Error 3%: 158 PO pairs, 111 other,s 80 trios. 47 SNPs.
# Error 5%: 133 PO pairs, 126 others, 85 trios. 47 SNPs. 

# conclusions here are as follows:
  # 1. even when using a highly tailored dataset, 47 SNPs is NOT enough for this
  # group of individuals at 5% error.

  # 2. Even just a few more quality SNPs allows for many more trios to be made.
  # (e.g., gmr_min8_thin5M_nomiss (47 SNPs) vs gmr5M_nomiss (43 SNPs) at the same
  # error rate (1%), those extra 4 SNPs (or some other attribute of those data)
  # allow for an extra 9 trios to be made.

  # 3. with high error (5%), you can overshoot the true number of trios by a big
  # margain using only ~50 SNPs.






##### Now for miss95 (missing data per site = 5% or less) #####
error_rate <- c(0.05, 0.05, 0.05)
gmr_min8_thin1M <- GetMaybeRel(GenoM = check_min8_thin1M,
                                      Err = error_rate,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin1M))
# Error 1%: 165 PO pairs, 17 others, 70 trios. 418 SNPs.
# Error 2.5%: 169 PO pairs, 21 others, 88 trios. 418 SNPs.
# Error 5%: 153 PO pairs, 45 others, 89 trios. 418 SNPs.

gmr_min8_thin2M <- GetMaybeRel(GenoM = check_min8_thin2M,
                               Err = error_rate,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin2M))
# Error 1%: 158 PO pairs, 24 others, 71 trios. 272 SNPs.
# Error 2.5%: 160 PO pairs, 35 others, 86 trios. 272 SNPs.
# Error 5%: 114 PO pairs, 61 others, 92 trios. 272 SNPs.

gmr_min8_thin3M <- GetMaybeRel(GenoM = check_min8_thin3M,
                               Err = error_rate,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin3M))
# Error 1%: 144 PO pairs, 36 others, 70 trios. 211 SNPs.
# Error 2.5%: 146 pairs, 40 others, 75 trios. 211 SNPs.
# Error 5%: 137 PO pairs, 85 others, 96 trios. 211 SNPs.

gmr_min8_thin4M <- GetMaybeRel(GenoM = check_min8_thin4M,
                               Err = error_rate,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin4M))
# Error 1%: 137 PO pairs, 46 others, 69 trios. 175 SNPs.
# Error 2.5%: 145 pairs, 51 others, 82 trios. 175 SNPs.
# Error 5%: 132 PO pairs, 92 others, 94 trios. 175 SNPs.

gmr_min8_thin5M <- GetMaybeRel(GenoM = check_min8_thin5M,
                               Err = error_rate,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin5M))
# Error 1%: 138 pairs, 43 others, 65 trios. 151 SNPs.
# Error 2.5%: 139 pairs, 53 others, 79 trios. 151 SNPs.
# Error 5%: 143 pairs, 100 others, 101 trios. 151 SNPs.

# This tells me that it's possible to get the right amount of trios with that many
# loci, it just requires increasing the error. I would imagine that there is some
# considerable error, and even with 5% error specified, we aren't getting all of 
# them because sequoia's too strict. 
# We should keep in mind, however, that this dataset DOES have missing data. Allowing
# no missing data, and some error rates with less than 100 SNPs allows us to come up
# with the 95 trios. HOWEVER, our ability to get sites with no missing data will decrease
# with more individuals added to the dataset. This might be something to reach out
# to those guys about. 
# ALL OF THAT SAID: It's still important to make sure that we're coming up with
# the CORRECT trios when we do hit that ~95 number. I think the most reasonable next
# step is to take a look at those gmr_mindep8_thinXM_nomiss objects and run the
# assignments thru the lookup table. This is to be done first thing 8/15.


#################################################################################

# Work 8/15/25 is going to start with trying a more refined error matrix. If that
# helps, I'm going to then start checking the assignments against the cross lookup
# table, but I'm not yet convinced that I'm there yet. 

erm1 <- ErrToM(Err = 0.01)
erm25 <- ErrToM(Err = 0.025)
erm5 <- ErrToM(Err = 0.05)
gmr_min8_thin1M_nomiss <- GetMaybeRel(GenoM = check_min8_thin1M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin1M_nomiss))
# erm1: 159 pairs, 31 others, 76 trios. 137 SNPs.
# erm25: 159 pairs, 34 others, 80 trios. 137 SNPs.
# erm5: 155 pairs, 49 others, 87 trios. 137 SNPs.

gmr_min8_thin2M_nomiss <- GetMaybeRel(GenoM = check_min8_thin2M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin2M_nomiss))
# erm1: 145 pairs, 50 others, 71 trios. 87 SNPs.
# erm25: 151 pairs, 61 others, 81 trios. 87 SNPs.
# erm5: 161 pairs, 85 others, 97 trios. 87 SNPs.

gmr_min8_thin3M_nomiss <- GetMaybeRel(GenoM = check_min8_thin3M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin3M_nomiss))
# erm1: 145 pairs, 50 others, 70 trios. 79 SNPs.
# erm25: 146 pairs, 68 others, 81 trios. 79 SNPs.
# erm5: 153 pairs, 115 others, 101 trios. 79 SNPs.

gmr_min8_thin4M_nomiss <- GetMaybeRel(GenoM = check_min8_thin4M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin4M_nomiss))
# erm1: 178 pairs, 75 others, 72 trios. 53 SNPs.
# erm25: 179 pairs, 98 others, 88 trios. 53 SNPs.
# erm5: 194 pairs, 131 others, 122 trios. 53 SNPs.

gmr_min8_thin5M_nomiss <- GetMaybeRel(GenoM = check_min8_thin5M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin5M_nomiss))
# erm1: 189 pairs, 68 others, 67 trios. 47 SNPs.
# erm25: 185 pairs, 94 others, 77 trios. 47 SNPs.
# erm5: 178 pairs, 134 others, 92 trios. 47 SNPs.

# with missing data!
gmr_min8_thin1M <- GetMaybeRel(GenoM = check_min8_thin1M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin1M))
# erm1: 141, 30, 57, 418.
# erm25: 156, 24, 71, 418.
# erm5: 163, 24, 85, 418.

gmr_min8_thin2M <- GetMaybeRel(GenoM = check_min8_thin2M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin2M))
# erm1: 131, 43, 59, 272.
# erm25: 148, 33, 74, 272.
# erm5: 152, 36, 81, 272.

gmr_min8_thin3M <- GetMaybeRel(GenoM = check_min8_thin3M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin3M))
# erm1: 118, 46, 54, 211.
# erm25: 134, 41, 66, 211.
# erm5: 138, 46, 71, 211.

gmr_min8_thin4M <- GetMaybeRel(GenoM = check_min8_thin4M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin4M))
# erm1: 117, 54, 55, 175.
# erm25: 129, 56, 72, 175.
# erm5: 133, 64, 82, 175.

gmr_min8_thin5M <- GetMaybeRel(GenoM = check_min8_thin5M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin5M))
# erm1: 109, 50, 52, 151.
# erm25: 124, 52, 65, 151.
# erm5: 130, 62, 76, 151.

################################################################################
# Looks like gmr_min8_thin2M_nomiss is our best one yet.
# # erm5: 161 pairs, 85 others, 97 trios. 87 SNPs.
table(unique(gmr_min8_thin2M_nomiss[["MaybeTrio"]]$id))
# missing 6703, 6709, 6711, 6722, 6725, 6736, 6747, 6756, 6757 (wasn't sequenced),
# 6764, 6782, 6788, 
# okay that sucks that i'm missing 11 individuals (offspring)

# now we need to see if the larger ones are even including all of the expected inds
# now we'll gmr_min8_thin3M_nomiss, erm5: 153 pairs, 115 others, 101 trios. 79 SNPs.
table(unique(gmr_min8_thin3M_nomiss[["MaybeTrio"]]$id))
# missing 6709, 6711, 6716, 6722, 6723, 6725, 6725, 6728, 6735, 6739, 6748, 6749,
# 6751, 6757 (not sequenced), 6764.
# missing 14 individuals on this one.

# that is absurd. i now need to see if the largest one i've come up with, which 
# is 122 trios (i think), captures all of the F1's. if not, we've got problems.
# if so, then it's just another thing to tweak and dial in from an LLR filtering
# or assignment threshold standpoint.
# gmr_min8_thin4M_nomiss 194 pairs, 131 others, 122 trios.
table(unique(gmr_min8_thin4M_nomiss[["MaybeTrio"]]$id))
# missing 6707, 6721, 6722, 6725, 6727, 6728, 6735, 6741, 6742, 6743, 6747, 6756, 
# 6757 (not sequenced), 6766, 6775, 6788
# MISSING 15 INDIVIDUALS HOW IS THAT EVEN POSSIBLE
table(unique(gmr_min8_thin4M_nomiss[["MaybePar"]]$ID1))
# missing 6735, 6757 (not sequenced)
# alright, so those individuals get parents placed to them, but just don't get trios.

################################################################################

# Work 8/22/25 -- Older Version 2.8.3

# because we can't get all of the test f1 individuals to be listed in the gmr objects
# we may have a problem with sequoia discriminating between two trios that have
# similar LLR or OH (see Jisca's email). she suggested using an older program version
# (e.g., 2.8.3), so that's what I'm trying today. first step is to get that older 
# version installed from github. 

# https://github.com/JiscaH/sequoia_archives/tree/main/2.8.3

install.packages("/Users/samjohnson/Desktop/sequoia_2.8.3_r-release_arm64.tgz", 
                 repos = NULL, type = "binary")
library(sequoia)

# these are the things you learn. don't be afraid to ask.
# now you've got it here.

erm1 <- ErrToM(Err = 0.01, flavour = "version2.0", Return = "matrix")
erm25 <- ErrToM(Err = 0.025, flavour = "version2.0", Return = "matrix")
erm5 <- ErrToM(Err = 0.05, flavour = "version2.0", Return = "matrix")

# error matrices established.
# next to do is to check the architecture of my gmr lines before use and to run 
# them on the error matrices above. also want to check the structure of the LH
# data against cran to make sure that i have the right kind of data going in there.

gmr_min8_thin1M_nomiss <- GetMaybeRel(GenoM = check_min8_thin1M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin1M_nomiss))
#erm1: 159 pairs, 31 others, 76 trios, 137 SNPs.
#erm25: 159 pairs, 34 others, 80 trios. 137 SNPs.
#erm5: 155 pairs, 49 others, 86 trios. 137 SNPs.

gmr_min8_thin2M_nomiss <- GetMaybeRel(GenoM = check_min8_thin2M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin2M_nomiss))
#erm1: 145 pairs, 50 others, 71 trios. 87 SNPs.
#erm25: 151 pairs, 61 others, 81 trios. 87 SNPs.
#erm5: 161 pairs, 85 others, 97 trios. 87 SNPs.

gmr_min8_thin3M_nomiss <- GetMaybeRel(GenoM = check_min8_thin3M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin3M_nomiss))
#erm1: 145 pairs, 50 others, 70 trios. 79 SNPs.
#erm25: 146 pairs, 68 others,  81 trios. 79 SNPs.
#erm5: 152 pairs, 115 others, 100 trios. 79 SNPs.

gmr_min8_thin4M_nomiss <- GetMaybeRel(GenoM = check_min8_thin4M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin4M_nomiss))
#erm1: 178 pairs, 75 others, 72 trios. 53 SNPs.
#erm25: 179 pairs, 98 others, 88 trios. 53 SNPs. 
#erm5: 194 pairs, 129 others, 121 trios. 53 SNPs.

gmr_min8_thin5M_nomiss <- GetMaybeRel(GenoM = check_min8_thin5M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin5M_nomiss))
#erm1: 189 pairs, 68 others, 67 trios. 47 SNPs.
#erm25: 185 pairs, 93 others, 77 trios. 47 SNPs.
#erm5: 177 pairs, 134 others, 91 trios. 47 SNPs. 

# with missing data!
gmr_min8_thin1M <- GetMaybeRel(GenoM = check_min8_thin1M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin1M))
#erm1: 141 pairs, 30 others, 57 trios. 418 SNPs.
#erm25: 156 pairs, 25 others, 72 trios. 418 SNPs.
#erm5: 163 pairs, 24 others, 85 trios. 418 SNPs.


gmr_min8_thin2M <- GetMaybeRel(GenoM = check_min8_thin2M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin2M))
#erm1: 131 pairs, 43 others, 59 trios. 272 SNPs.
#erm25: 148 pairs, 33 others, 73 trios. 272 SNPs.
#erm5: 152 pairs, 36 others, 80 trios. 272 SNPs.

gmr_min8_thin3M <- GetMaybeRel(GenoM = check_min8_thin3M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin3M))
#erm1: 118 pairs, 46 others, 54 trios. 211 SNPs.
#erm25: 134 pairs, 41 others, 66 trios. 211 SNPs.
#erm5: 138 pairs, 47 others, 72 trios. 211 SNPs. 


gmr_min8_thin4M <- GetMaybeRel(GenoM = check_min8_thin4M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin4M))
#erm1: 117 pairs, 54 others, 55 trios. 175 SNPs.
#erm25: 129 pairs, 56 others, 72 trios. 175 SNPs.
#erm5: 133 pairs, 64 others, 82 trios. 175 SNPs. 


gmr_min8_thin5M <- GetMaybeRel(GenoM = check_min8_thin5M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin5M))
#erm1: 109 pairs, 50 others, 52 trios. 151 SNPs.
#erm25: 124 pairs, 52 others, 65 trios. 151 SNPs.
#erm5: 129 pairs, 62 others, 75 trios. 151 SNPs. 


# these three are the only ones that surpassed 95 trios. 
gmr_min8_thin2M_nomiss #97
gmr_min8_thin3M_nomiss #100
gmr_min8_thin4M_nomiss #121

# and on the old version, it's been 97, 101, and 122. so we're LOSING in some cases.
# maybe they're better assignments. let's look at composition.
table(unique(gmr_min8_thin2M_nomiss[["MaybeTrio"]]$id))
length(unique(gmr_min8_thin2M_nomiss[["MaybeTrio"]]$id))
# 82 indivs
length(unique(gmr_min8_thin2M_nomiss[["MaybePar"]]$ID1))
# 95 individuals listed in the MaybePar, just not trios. 

table(unique(gmr_min8_thin3M_nomiss[["MaybeTrio"]]$id))
length(unique(gmr_min8_thin3M_nomiss[["MaybeTrio"]]$id))
# 81 indivs
length(unique(gmr_min8_thin3M_nomiss[["MaybePar"]]$ID1))
# 94 individuals listed in the MaybePar, just not trios.

table(unique(gmr_min8_thin4M_nomiss[["MaybeTrio"]]$id))
length(unique(gmr_min8_thin4M_nomiss[["MaybeTrio"]]$id))
# 80 indivs
length(unique(gmr_min8_thin4M_nomiss[["MaybePar"]]$ID1))
# 94 individuals listed in the MaybePar, just not trios.

# what is happening to the others of those tests?

# next steps. get trios we feel okay about, look at dist of LLR and OH for
# those and then expand? just don't know what else to do.

# wonder if the problem has to do with the fact that i'm operating on Complex = "simp"
# which says that no explicit consideration of inbred relationships

gmr_min8_thin2M_nomiss <- GetMaybeRel(GenoM = check_min8_thin2M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "full",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin2M_nomiss))
#complex = "simp"
#erm1: 145 pairs, 50 others, 71 trios. 87 SNPs.
#erm25: 151 pairs, 61 others, 81 trios. 87 SNPs.
#erm5: 161 pairs, 85 others, 97 trios. 87 SNPs.

#complex = "full"
#erm1: 144 pairs, 42 others, 70 trios. 87 SNPs.
#erm25: 148 pairs, 50 others, 78 trios. 87 SNPs.
#erm5: 157 pairs, 67 others, 90 trios. 87 SNPs.


gmr_min8_thin3M_nomiss <- GetMaybeRel(GenoM = check_min8_thin3M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "full",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin3M_nomiss))
#complex = "simp"
#erm1: 145 pairs, 50 others, 70 trios. 79 SNPs.
#erm25: 146 pairs, 68 others,  81 trios. 79 SNPs.
#erm5: 152 pairs, 115 others, 100 trios. 79 SNPs.

#complex = "full"
#erm1: 142 pairs, 43 others, 68 trios. 79 SNPs.
#erm25: 144 pairs, 57 others, 78 trios. 79 SNPs.
#erm5: 150 pairs, 93 others, 95 trios. 79 SNPs.

gmr_min8_thin4M_nomiss <- GetMaybeRel(GenoM = check_min8_thin4M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "full",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin4M_nomiss))
#complex = "simp"
#erm1: 178 pairs, 75 others, 72 trios. 53 SNPs.
#erm25: 179 pairs, 98 others, 88 trios. 53 SNPs. 
#erm5: 194 pairs, 129 others, 121 trios. 53 SNPs.

#complex = "full"
#erm1: 179 pairs, 69 others, 73 trios. 53 SNPs.
#erm25: 178 pairs, 84 others, 87 trios. 53 SNPs.
#erm5: 191 pairs, 113 others, 116 trios. 53 SNPs.

table(unique(gmr_min8_thin3M_nomiss[["MaybeTrio"]]$id))
length(unique(gmr_min8_thin3M_nomiss[["MaybeTrio"]]$id))
# 79 indivs
length(unique(gmr_min8_thin3M_nomiss[["MaybePar"]]$ID1))
# 94 individuals listed in the MaybePar, just not trios.

table(unique(gmr_min8_thin4M_nomiss[["MaybeTrio"]]$id))
length(unique(gmr_min8_thin4M_nomiss[["MaybeTrio"]]$id))
# 79 indivs
length(unique(gmr_min8_thin4M_nomiss[["MaybePar"]]$ID1))
# 94 individuals listed in the MaybePar, just not trios.

# difference is negligable. wonder what's going on here. going to try more rigourous
# LD filtering to start i think. genotype filtering. check out the ageprior. per snp 
# error rate?

ped_check_min8_thin2M_nomiss <- sequoia(GenoM = check_min8_thin2M_nomiss,
                                        LifeHistData =  LH_Data,
                                        Module = "par",
                                        Err = erm5,
                                        #Tfilter =, 
                                        Tassign = 1.0,
                                        Complex = "simp",
                                        UseAge = "yes")
gmr_ped_check_min8_thin2M_nomiss <- GetMaybeRel(GenoM = check_min8_thin4M_nomiss,
                                                Err = erm5,
                                                SeqList = ped_check_min8_thin2M_nomiss,
                                                Module = "par",
                                                # MaxMismatch = NA,
                                                Complex = "full",
                                                LifeHistData = LH_Data,
                                                quiet = FALSE,
                                                Tassign = 1.0,
                                                # Tfilter = -100,
                                                MaxPairs = 7*nrow(check_min8_thin4M_nomiss))
# 82 pairs, 0 others, 20 trios. what.

# i'm also getting confused as to what's up with the assignment threshold.
# alright, check this out. when i take Tassign to 0.01

gmr_min8_thin4M_nomiss <- GetMaybeRel(GenoM = check_min8_thin4M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "full",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 0.01,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin4M_nomiss))

# 270 pairs, 220 others, 137 trios. 53 SNPs. 
length(unique(gmr_min8_thin4M_nomiss[["MaybePar"]]$ID1))
table(unique(gmr_min8_thin4M_nomiss[["MaybePar"]]$ID1))
length(unique(gmr_min8_thin4M_nomiss[["MaybeTrio"]]$id))
table(unique(gmr_min8_thin4M_nomiss[["MaybeTrio"]]$id))
# alright, up to 84 unique individuals (offspring) at the center of those trios.

# i've noticed there that in MaybePar, some of the FS relationships are between
# test F1's and F0's. impossible. does that get weeded out when we add a min and
# max birth year to the LH_Data?
LH_Data <- data.frame(LH_Data, BY.min = NA, BY.max = NA)
LH_Data[1:95, 4:5] <- 2005
LH_Data[96:208, 4:5] <- 2000
#LH_Data <- LH_Data %>% 
    #select(!BirthYear)

gmr_min8_thin4M_nomiss <- GetMaybeRel(GenoM = check_min8_thin4M_nomiss,
                                      Err = erm5,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "full",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin4M_nomiss))
# 191, 113 others, 116 trios.
length(unique(gmr_min8_thin4M_nomiss[["MaybeTrio"]]$id))
# 79 unique indivs.

# still doing the thing in the MaybePar where FS have an AgeDif of 5, so let's
# go ahead and explore MakeAgePrior

ped_check_min8_thin2M_nomiss <- sequoia(GenoM = check_min8_thin4M_nomiss,
                                        LifeHistData =  LH_Data,
                                        Module = "par",
                                        Err = erm5,
                                        #Tfilter =, 
                                        Tassign = 1.0,
                                        Complex = "simp",
                                        UseAge = "yes")



################################################################################

# Work 8/25/25 -- Newest Version 3.0.3 : Just released! Whoop:

# Jisca says "Please have a try with the newest version I released on CRAN earlier 
# this week, it has a whole bunch of big fixes and I tuned it to hopefully have 
# a better balance between false positives and false negatives."

install.packages("sequoia")
library(sequoia)

erm1 <- ErrToM(Err = 0.01, flavour = "version2.0", Return = "matrix")
erm25 <- ErrToM(Err = 0.025, flavour = "version2.0", Return = "matrix")
erm5 <- ErrToM(Err = 0.05, flavour = "version2.0", Return = "matrix")

# i'm still concerned with my LH_Data, and the absence of min/max birth years.
# let's try it with what we've got now, and see what happens. plenty of other
# upstream options to investigate.

gmr_min8_thin1M_nomiss <- GetMaybeRel(GenoM = check_min8_thin1M_nomiss,
                                      Err = erm1,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "full",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 0.01,
                                      Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin1M_nomiss))
#erm1: 159, 31, 95, 137.
length(unique(gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id))
table(unique(gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id))
# 85 unique indivs

#erm25: 159, 34, 127, 137.
length(unique(gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id))
# 91 unique indivs!!!!

#erm5: 155, 49, 184, 137. 
length(unique(gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id))
# 109 unique indivs... alright, something's fishy now.

# it's assigning some of the F0 individuals to F1 parents... ~_~ why...
# also still assigning F0 F1 FS pairs.

# let's alter the BY's and see if that helps.
LH_Data <- data.frame(LH_Data, BY.min = NA, BY.max = NA)
LH_Data[1:95, 4:5] <- 2001
LH_Data[96:208, 4:5] <- 2000
LH_Data[1:95, 3] <- 2001
LH_Data[96:208, 3] <- 2000
LH_Data$BirthYear <- NA
LH_Data$Sex[LH_Data$Sex == "M"] <- 2
LH_Data$Sex[LH_Data$Sex == "F"] <- 1
LH_Data$Sex <- as.numeric(LH_Data$Sex)
#LH_Data <- LH_Data %>% 
  #select(!BirthYear)
# BY.min and BY.max are ignored if BirthYear exists.

# see output from gmr, default MaxAgeParent = 99,99
# set LH data accordingly ^^^
# now erm1: gmr_min8_thin1M_nomiss is 189, 0, 101 trios, 137.
length(unique(gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id))
# 89 indivs
table(unique(gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id))
# still with some F0's as the offspring.

# --- All of this was unsuccessfull work 08/26/25, I just don't yet want to get
#     rid of it in case it proves valuable later --- #
# age_prior <- MakeAgePrior(LifeHistData = LH_Data, 
#                           MinAgeParent = 1,
#                           MaxAgeParent = 1,
#                           Discrete = 1)
# 
# gmr_min8_thin1M_nomiss <- GetMaybeRel(GenoM = check_min8_thin1M_nomiss,
#                                       Err = erm1,
#                                       # SeqList = outfull,
#                                       Module = "par",
#                                       # MaxMismatch = NA,
#                                       Complex = "simp",
#                                       LifeHistData = LH_Data,
#                                       AgePrior = age_prior,
#                                       quiet = FALSE,
#                                       Tassign = 1,
#                                       #Tfilter = -10,
#                                       MaxPairs = 7*nrow(check_min8_thin1M_nomiss))
# length(unique(gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id))
# table(unique(gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id))
# 
# 
# 
# maf <- sequoia::SnpStats(check_min8_thin1M_nomiss, Plot=FALSE)[,"AF"]
# erm2 <- ErrToM(Err = 0.02, flavour="version2.9", Return="matrix")
# mm <- CalcMaxMismatch(Err = 0.02, MAF = maf)  # returns DUP, OH, ME
# mm
# 
# mm["OH"] <- mm["OH"] + 3
# mm["ME"] <- mm["ME"] + 5
# 
# 
# ped_min8_thin1M_nomiss <- sequoia(
#   GenoM   = check_min8_thin1M_nomiss,
#   Err     = erm2,
#   Module  = "par",
#   Complex = "simp",
#   LifeHistData = LH_Data,
#   UseAge  = "yes",
#   args.AP = list(
#     MinAgeParent = 1,
#     MaxAgeParent = 1,
#     Discrete = 1),
#   quiet   = FALSE,
#   Tassign = 0.5,
#   Tfilter = -10,
#   StrictGenoCheck = TRUE)

##### Work 08/26/25 - Checking all of my existing matrices in version 3.0.3 #####
## objective is to report the number of unique test F1's that are in the trios ##

# LET'S GET CRAZY
erm75 <- ErrToM(Err = 0.075)
erm10 <- ErrToM(Err = 0.1)

gmr_min8_thin1M_nomiss <- GetMaybeRel(GenoM = check_min8_thin1M_nomiss,
                                      Err = erm1,
                                      #SeqList = ped_min8_thin1M_nomiss,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      #AgePrior = age_prior,
                                      quiet = FALSE,
                                      Tassign = 1,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin1M_nomiss))
# no age prior
#erm1: 159, 31, 95, 137. 76 unique f1's in the trios.
#erm25: 159, 34, 127, 137. 79 test f1's in the trios.
#erm5: 155, 49, 184, 137. 83 test f1's in the trios.
#erm75: 146, 75, 208, 137. 71 test f1's in the trios.
#erm10: 146, 114, 208, 137. 27 test f1's in the trios.

# here's how you tell which f1's, and how many
unique(gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id)])
length(unique(gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id)]))

# how many unique focal indivs are in the trios
length(unique(gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id))
# see who it is
table(unique(gmr_min8_thin1M_nomiss[["MaybeTrio"]]$id))



gmr_min8_thin2M_nomiss <- GetMaybeRel(GenoM = check_min8_thin2M_nomiss,
                                      Err = erm10,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin2M_nomiss))
# no age prior
#erm1: 145, 50, 90, 87. 69 test f1's in the trios.
#erm25: 151, 61, 131. 87. 78 test f1's in the trios.
#erm5: 161, 85, 200. 87. 82 test f1's in the trios.
#erm75: 164, 128, 208, 87. 48 test f1's in the trios.
#erm10: 162, 179, 208, 87. 3 test f1's in the trios.

length(unique(gmr_min8_thin2M_nomiss[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_min8_thin2M_nomiss[["MaybeTrio"]]$id)]))



gmr_min8_thin3M_nomiss <- GetMaybeRel(GenoM = check_min8_thin3M_nomiss,
                                      Err = erm10,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin3M_nomiss))
# no age prior
#erm1: 145, 50, 88, 79. 67 test f1's in the trios.
#erm25: 146, 68, 129, 79. 75 test f1's in the trios.
#erm5: 152, 115, 180, 79. 81 test f1's in the trios.
#erm75: 152, 175, 208, 79. 53 test f1's in the trios.
#erm10: 157, 239, 208, 79. 9 test f1's in the trios. 

length(unique(gmr_min8_thin3M_nomiss[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_min8_thin3M_nomiss[["MaybeTrio"]]$id)]))



gmr_min8_thin4M_nomiss <- GetMaybeRel(GenoM = check_min8_thin4M_nomiss,
                                      Err = erm75,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin4M_nomiss))
#erm1: 178, 75, 114, 53. 62 test f1's in the trios.
#erm25: 179, 98, 157, 53. 69 test f1's in the trios.
#erm5: 194, 129, 208, 53. 54 test f1's in the trios.
#erm75: 195, 183, 208, 53. 4 test f1's in the trios. 
#erm10: 188, 228, 208, 53. 0 test f1's in the trios.

length(unique(gmr_min8_thin4M_nomiss[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_min8_thin4M_nomiss[["MaybeTrio"]]$id)]))



gmr_min8_thin5M_nomiss <- GetMaybeRel(GenoM = check_min8_thin5M_nomiss,
                                      Err = erm10,
                                      # SeqList = outfull,
                                      Module = "par",
                                      # MaxMismatch = NA,
                                      Complex = "simp",
                                      LifeHistData = LH_Data,
                                      quiet = FALSE,
                                      Tassign = 1.0,
                                      # Tfilter = -100,
                                      MaxPairs = 7*nrow(check_min8_thin5M_nomiss))
#erm1: 189, 68, 107, 47. 58 test f1's in the trios.
#erm25: 185, 93, 138, 47. 65 test f1's in the trios.
#erm5: 177, 134, 204, 47. 67 test f1's in the trios.
#erm75: 168, 181, 208, 47. 28 test f1's in the trios.
#erm10: 165, 233, 208, 47. 5 test f1's in the trios.

length(unique(gmr_min8_thin5M_nomiss[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_min8_thin5M_nomiss[["MaybeTrio"]]$id)]))

# alright, maybe with missing data and more SNPs?!

################################################################################

# with missing data!
gmr_min8_thin1M <- GetMaybeRel(GenoM = check_min8_thin1M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin1M))
#erm1: 141, 30, 60, 418. 57 test f1's in the trios.
#erm25: 156, 25, 98, 418. 72 test f1's in the trios.
#erm5: 163, 24, 160, 418. 85 test f1's in the trios. 

length(unique(gmr_min8_thin1M[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_min8_thin1M[["MaybeTrio"]]$id)]))



gmr_min8_thin2M <- GetMaybeRel(GenoM = check_min8_thin2M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin2M))
#erm1: 131, 43, 68, 272. 59 test f1's in the trios.
#erm25: 148, 33, 106, 272. 73 test f1's in the trios.
#erm5: 152, 36, 159, 272. 80 test f1's in the trios.

length(unique(gmr_min8_thin2M[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_min8_thin2M[["MaybeTrio"]]$id)]))



gmr_min8_thin3M <- GetMaybeRel(GenoM = check_min8_thin3M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin3M))
#erm1: 118, 46, 66, 211. 54 test f1's in the trios.
#erm25: 134, 41, 95, 211. 66 test f1's in the trios.
#erm5: 138, 47, 139, 211. 72 test f1's in the trios.

length(unique(gmr_min8_thin3M[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_min8_thin3M[["MaybeTrio"]]$id)]))



gmr_min8_thin4M <- GetMaybeRel(GenoM = check_min8_thin4M,
                               Err = erm5,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin4M))
#erm1: 117, 54, 62, 175. 55 test f1's in the trios.
#erm25: 129, 56, 92, 175. 72 test f1's in the trios.
#erm5: 133, 64, 140, 175. 81 test f1's in the trios.

length(unique(gmr_min8_thin4M[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_min8_thin4M[["MaybeTrio"]]$id)]))



gmr_min8_thin5M <- GetMaybeRel(GenoM = check_min8_thin5M,
                               Err = erm25,
                               # SeqList = outfull,
                               Module = "par",
                               # MaxMismatch = NA,
                               Complex = "simp",
                               LifeHistData = LH_Data,
                               quiet = FALSE,
                               Tassign = 1.0,
                               # Tfilter = -100,
                               MaxPairs = 7*nrow(check_min8_thin5M))
#erm1: 109, 50, 64, 151. 52 test f1's in the trios.
#erm25: 124, 52, 92. 151. 65 test f1's in the trios.
#erm5: 129, 62, 139, 151. 75 test f1's in the trios.

length(unique(gmr_min8_thin5M[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr_min8_thin5M[["MaybeTrio"]]$id)]))

#####  ----- Alright, so we got to 85/95 with gmr_min8_thin1M with erm5 ----- #####

trios_gmr_min8_thin1M_erm5 <- gmr_min8_thin1M[["MaybeTrio"]]
head(trios_gmr_min8_thin1M_erm5)

# read in lookup table
cross_lookup <- read.csv(file = "true_cross_lookup.csv", header = TRUE)
cross_lookup <- cross_lookup[,-1]

# this is a ridiculously effective set of piped functions here:
trios_gmr_min8_thin1M_erm5 <- trios_gmr_min8_thin1M_erm5 %>%
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

head(trios_gmr_min8_thin1M_erm5)
table(trios_gmr_min8_thin1M_erm5$valid_cross)


trios_gmr_min8_thin1M_erm5_checked <- trios_gmr_min8_thin1M_erm5 %>%
  select(id, parent1, parent2, pair, valid_cross, everything())

trios_gmr_min8_thin1M_erm5_checked <- trios_gmr_min8_thin1M_erm5_checked %>% 
  filter(!LLRparent1 %in% c(555, -555),
         !LLRparent2 %in% c(555, -555),
         !LLRpair    %in% c(555, -555))
dim(trios_gmr_min8_thin1M_erm5_checked)
table(trios_gmr_min8_thin1M_erm5_checked$valid_cross) # 76 true, 9 false. 89.41% 

#####  ----- And we got to 81/95 with gmr_min8_thin4M with erm5 ----- #####

trios_gmr_min8_thin4M_erm5 <- gmr_min8_thin4M[["MaybeTrio"]]
head(trios_gmr_min8_thin4M_erm5)

# this is a ridiculously effective set of piped functions here:
trios_gmr_min8_thin4M_erm5 <- trios_gmr_min8_thin4M_erm5 %>%
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

head(trios_gmr_min8_thin4M_erm5)
table(trios_gmr_min8_thin4M_erm5$valid_cross)


trios_gmr_min8_thin4M_erm5_checked <- trios_gmr_min8_thin4M_erm5 %>%
  select(id, parent1, parent2, pair, valid_cross, everything())

trios_gmr_min8_thin4M_erm5_checked <- trios_gmr_min8_thin4M_erm5_checked %>% 
  filter(!LLRparent1 %in% c(555, -555),
         !LLRparent2 %in% c(555, -555),
         !LLRpair    %in% c(555, -555))
dim(trios_gmr_min8_thin4M_erm5_checked)
table(trios_gmr_min8_thin4M_erm5_checked$valid_cross) # 73 true, 9 false. 89.02% 
# 175 SNPs, up to 12 ME mismatches. That's 6.8%. What if we crank the error?

# go check the no missing data. error rates over 5% are a no go (as if that wasn't
# already obvious from the documentation)

# Seems like the success rate is relatively similar to hiphop, which is good.
# Somewhere in the ballpack of 90%, it would just be nice to know that I'm picking
# up more crosses than what I'm getting right now. Figure that adding SNPs might
# be a good thing, as well as more rigorous filtering to LD and the genotypes themselves.

##### ----- 08/27/25 Work: r^2 coefficients (LD) for within-scaffold pairs ----- #####

rsq_500kb <- read.table(file = "nothin_window500K.geno.ld", header = TRUE)
rsq_100kb <- read.table(file = "nothin_window100K.geno.ld", header = TRUE)

rsq_500kb$dist <- abs(rsq_500kb$POS2 - rsq_500kb$POS1)
rsq_100kb$dist <- abs(rsq_100kb$POS2 - rsq_100kb$POS1)

rsq_500kb$CHR <- as.factor(rsq_500kb$CHR)
rsq_100kb$CHR <- as.factor(rsq_100kb$CHR)

rsq_500kb$CHR_num <- as.numeric(sub("scaffold_","", rsq_500kb$CHR))
rsq_100kb$CHR_num <- as.numeric(sub("scaffold_","", rsq_100kb$CHR))

cor.test(rsq_500kb$dist, rsq_500kb$R.2, method = "pearson")
cor.test(rsq_500kb$dist, rsq_500kb$R.2, method = "spearman")

cor.test(rsq_100kb$dist, rsq_100kb$R.2, method = "pearson")
cor.test(rsq_100kb$dist, rsq_100kb$R.2, method = "spearman")


# Plot 500Kb
#------------------------------------------------------------------------------#
scaf_levels_500 <- sort(unique(rsq_500kb$CHR_num))
n_scaf_500 <- length(scaf_levels_500)
mycols_500 <- rainbow(n_scaf_500)
rsq_500kb$col <- mycols_500[match(rsq_500kb$CHR_num, scaf_levels_500)]

par(mar = c(5, 4, 4, 8)) 
plot(R.2 ~ dist, data = rsq_500kb, col = col, pch = 19, 
     main = "LD by distance (within-scaffold snps, 500kb window)")

legend("topright", 
       inset = c(-0.15, 0),
       legend = paste0("scaffold_", scaf_levels),
       col = mycols,
       pch = 19,
       cex = 0.7,
       xpd = TRUE)
#------------------------------------------------------------------------------#




# Plot 100Kb
#------------------------------------------------------------------------------#
scaf_levels_100 <- sort(unique(rsq_100kb$CHR_num))
n_scaf_100 <- length(scaf_levels_100)
mycols_100 <- rainbow(n_scaf_100)
rsq_100kb$col <- mycols_100[match(rsq_100kb$CHR_num, scaf_levels_100)]

par(mar = c(5, 4, 4, 8)) 
plot(R.2 ~ dist, data = rsq_100kb, col = col, pch = 19, 
     main = "LD by distance (within-scaffold snps, 100kb window)")

legend("topright", 
       inset = c(-0.15, 0),
       legend = paste0("scaffold_", scaf_levels),
       col = mycols,
       pch = 19,
       cex = 0.7,
       xpd = TRUE)
#------------------------------------------------------------------------------#



