## Sequoia_ContamFilt.R by SPJ 081125
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
# q20 - keep sites with quality > 20 (mean base quality? or is this mean q for the read?)
# mindep4 - minimum mean read depth 4 or above (removes low quality sites, low confidence calls)
# maxdep75 - maximum mean read depth 75 or lower (removes paralogs)
# maf30 - minor allele frequency 30% or higher
# miss95 - each snp needs to have 5% or less missing data

# here i will continue to implement the following list of suggestions from
# Dr. Jisca Huisman. In order of implementation, they are:
# (numbered by the scheme from her email)
# 3. filter the vcfs for LD by thinnning
      # all input matrices in this script are thinned by at least 1Mb
# 1. Increase the Tassign threshold
      
# 2. Alter the complex argument from sequoia()
  
# 4. Try an older sequoia program version (currently running 2.11.2)

# here's today's approach (08/11/25)
# read in these matrices. they're all thinned (1Mb, 2, 3, 4, 5)
# missing data has been an issue in the past, and sequoia says to limit md
# so let's remove sites with ANY missing data and see how many are remaining.

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

################################################################################
#### RUNNING ON FURTHER THINNED (LD FILTERED) FILES
################################################################################

setwd("/Users/samjohnson/Documents/Sauger_042225/GeneticData/Sequoia/Sequoia_Inp/contam_fastp_svit_mem/firstfilt_hardfilt_thinned/maf30/geno_mat")

# read in genotype matrices
mat1M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_mindep4_maxdep75_maf30_miss95_thin1M.012_conv", 
                    header = FALSE, sep = "\t", na.strings = c("NA", "-1"))
mat2M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_mindep4_maxdep75_maf30_miss95_thin2M.012_conv", 
                      header = FALSE, sep = "\t", na.strings = c("NA", "-1"))
mat3M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_mindep4_maxdep75_maf30_miss95_thin3M.012_conv", 
                    header = FALSE, sep = "\t", na.strings = c("NA", "-1"))
mat4M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_mindep4_maxdep75_maf30_miss95_thin4M.012_conv", 
                    header = FALSE, sep = "\t", na.strings = c("NA", "-1"))
mat5M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_mindep4_maxdep75_maf30_miss95_thin5M.012_conv", 
                    header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

# correct row names as sample id's (someday we should parallelize this)
mat1M <- mat1M[, -1]
mat2M <- mat2M[, -1]
mat3M <- mat3M[, -1]
mat4M <- mat4M[, -1]
mat5M <- mat5M[, -1]

mat1M <- as.matrix(mat1M)
mat2M <- as.matrix(mat2M)
mat3M <- as.matrix(mat3M)
mat4M <- as.matrix(mat4M)
mat5M <- as.matrix(mat5M)

# yes, i am aware that these are all the same. trying to be thorough. would hate 
# to mess up a matrix because i was lazy on a step like this.
ind1M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_mindep4_maxdep75_maf30_miss95_thin1M.012.indv", header = FALSE)
ind2M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_mindep4_maxdep75_maf30_miss95_thin2M.012.indv", header = FALSE)
ind3M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_mindep4_maxdep75_maf30_miss95_thin3M.012.indv", header = FALSE)
ind4M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_mindep4_maxdep75_maf30_miss95_thin4M.012.indv", header = FALSE)
ind5M <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q20_mindep4_maxdep75_maf30_miss95_thin5M.012.indv", header = FALSE)

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

rownames(mat1M) <- ind1M$sample
rownames(mat2M) <- ind2M$sample
rownames(mat3M) <- ind3M$sample
rownames(mat4M) <- ind4M$sample
rownames(mat5M) <- ind5M$sample

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
mat1M <- mat1M[rownames(mat1M) %in% LH_Data$ID, , drop = FALSE]
dim(mat1M) # 208 indivs, 437 loci

mat2M <- mat2M[rownames(mat2M) %in% LH_Data$ID, , drop = FALSE]
dim(mat2M) # 208 indivs, 280 loci

mat3M <- mat3M[rownames(mat3M) %in% LH_Data$ID, , drop = FALSE]
dim(mat3M) # 208 indivs, 215 loci

mat4M <- mat4M[rownames(mat4M) %in% LH_Data$ID, , drop = FALSE]
dim(mat4M) # 208 indivs, 182 loci

mat5M <- mat5M[rownames(mat5M) %in% LH_Data$ID, , drop = FALSE]
dim(mat5M) # 208 indivs, 154 loci

# check for all heterozygous sites since those likely will not be informative.
all_ones_1M <- apply(mat1M, 2, function(col) all(col == 1))
all_ones_2M <- apply(mat2M, 2, function(col) all(col == 1))
all_ones_3M <- apply(mat3M, 2, function(col) all(col == 1))
all_ones_4M <- apply(mat4M, 2, function(col) all(col == 1))
all_ones_5M <- apply(mat5M, 2, function(col) all(col == 1))

mat1M <- mat1M[, !all_ones_1M]
dim(mat1M) # none lost

mat2M <- mat2M[, !all_ones_2M]
dim(mat2M) # ONE lost

mat3M <- mat3M[, !all_ones_3M]
dim(mat3M) # none lost

mat4M <- mat4M[, !all_ones_4M]
dim(mat4M) # none lost

mat5M <- mat5M[, !all_ones_5M]
dim(mat5M) # none lost

# alright, looks like heterozygous sites are not going to be a problem here.
rownames(mat1M) <- ind1M$sample
rownames(mat2M) <- ind2M$sample
rownames(mat3M) <- ind3M$sample
rownames(mat4M) <- ind4M$sample
rownames(mat5M) <- ind5M$sample
# going back to this the next day, why does this chunk appear again down here?
# each of these matrices mat1-5M should already have their rownames filtered to
# include only the F0's and Test F1's...? Wonder why this is here.

# check genotype matrix for samples/loci to be excluded
check1M <- CheckGeno(mat1M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check2M <- CheckGeno(mat2M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check3M <- CheckGeno(mat3M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check4M <- CheckGeno(mat4M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
check5M <- CheckGeno(mat5M, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))

# alright, let's now remove any sites with missing data and see how many are retained
# create an index of sites where there are no -9 values
completesites_1M <- apply(check1M, 2, function(col) all(col != -9))
completesites_2M <- apply(check2M, 2, function(col) all(col != -9))
completesites_3M <- apply(check3M, 2, function(col) all(col != -9))
completesites_4M <- apply(check4M, 2, function(col) all(col != -9))
completesites_5M <- apply(check5M, 2, function(col) all(col != -9))

# subset the matrix to keep only those columns
check1M_nomiss <- check1M[, completesites_1M]
dim(check1M_nomiss) # 129 SNPs

check2M_nomiss <- check2M[, completesites_2M]
dim(check2M_nomiss) # 79 SNPs

check3M_nomiss <- check3M[, completesites_3M]
dim(check3M_nomiss) # 65 SNPs

check4M_nomiss <- check4M[, completesites_4M]
dim(check4M_nomiss) # 51 SNPs

check5M_nomiss <- check5M[, completesites_5M]
dim(check5M_nomiss) # 43 SNPs

# alright, now there are a few axes along which we can vary these datasets
# to analyze their usefulness. these are:
# number of loci, md, risk of complication due to LD (thinness), and error.
# the first two are intrinsic properties of the data themselves, while the third
# is a property of how sequoia will act on the data.

# first let's look at thinning. we know we can get ALMOST all of the pairs with
# as few as 64 SNPs, potentially fewer. so let's see how each of these matrices
# compare with the same error. 

error_rate <- c(0.01, 0.01, 0.01)
gmr1M_nomiss <- GetMaybeRel(GenoM = check1M_nomiss,
                            Err = error_rate,
                            # SeqList = outfull,
                            Module = "par",
                            # MaxMismatch = NA,
                            Complex = "simp",
                            LifeHistData = LH_Data,
                            quiet = FALSE,
                            Tassign = 1.0,
                            # Tfilter = -100,
                            MaxPairs = 7*nrow(check1M_nomiss))
# 157 PO pairs, 34 others, 80 trios. 129 SNPs.

gmr2M_nomiss <- GetMaybeRel(GenoM = check2M_nomiss,
                            Err = error_rate,
                            # SeqList = outfull,
                            Module = "par",
                            # MaxMismatch = NA,
                            Complex = "simp",
                            LifeHistData = LH_Data,
                            quiet = FALSE,
                            Tassign = 1.0,
                            # Tfilter = -100,
                            MaxPairs = 7*nrow(check2M_nomiss))
# 166 PO pairs, 58 others, 87 trios. 79 SNPs

gmr3M_nomiss <- GetMaybeRel(GenoM = check3M_nomiss,
                            Err = error_rate,
                            # SeqList = outfull,
                            Module = "par",
                            # MaxMismatch = NA,
                            Complex = "simp",
                            LifeHistData = LH_Data,
                            quiet = FALSE,
                            Tassign = 1.0,
                            # Tfilter = -100,
                            MaxPairs = 7*nrow(check3M_nomiss))
# 166 PO pairs, 77 others, 85 trios. 65 SNPs

gmr4M_nomiss <- GetMaybeRel(GenoM = check4M_nomiss,
                            Err = error_rate,
                            # SeqList = outfull,
                            Module = "par",
                            # MaxMismatch = NA,
                            Complex = "simp",
                            LifeHistData = LH_Data,
                            quiet = FALSE,
                            Tassign = 1.0,
                            # Tfilter = -100,
                            MaxPairs = 7*nrow(check4M_nomiss))
# 166 PO pairs, 92 others, 80 trios. 51 SNPs.

gmr5M_nomiss <- GetMaybeRel(GenoM = check5M_nomiss,
                            Err = error_rate,
                            # SeqList = outfull,
                            Module = "par",
                            # MaxMismatch = NA,
                            Complex = "simp",
                            LifeHistData = LH_Data,
                            quiet = FALSE,
                            Tassign = 1.0,
                            # Tfilter = -100,
                            MaxPairs = 7*nrow(check5M_nomiss))
# 170 PO pairs, 100 others, 63 trios. 43 SNPs. 

##### how does this compare when there is missing data? (miss95) #####
error_rate <- c(0.01, 0.01, 0.01)
gmr1M <- GetMaybeRel(GenoM = check1M,
                            Err = error_rate,
                            # SeqList = outfull,
                            Module = "par",
                            # MaxMismatch = NA,
                            Complex = "simp",
                            LifeHistData = LH_Data,
                            quiet = FALSE,
                            Tassign = 1.0,
                            # Tfilter = -100,
                            MaxPairs = 7*nrow(check1M))
# 158 PO pairs, 22 others, 64 trios. 437 SNPs.

gmr2M <- GetMaybeRel(GenoM = check2M,
                     Err = error_rate,
                     # SeqList = outfull,
                     Module = "par",
                     # MaxMismatch = NA,
                     Complex = "simp",
                     LifeHistData = LH_Data,
                     quiet = FALSE,
                     Tassign = 1.0,
                     # Tfilter = -100,
                     MaxPairs = 7*nrow(check2M))
# 147 PO pairs, 29 others, 65 trios. 279 SNPs.

gmr3M <- GetMaybeRel(GenoM = check3M,
                     Err = error_rate,
                     # SeqList = outfull,
                     Module = "par",
                     # MaxMismatch = NA,
                     Complex = "simp",
                     LifeHistData = LH_Data,
                     quiet = FALSE,
                     Tassign = 1.0,
                     # Tfilter = -100,
                     MaxPairs = 7*nrow(check3M))
# 133 PO pairs, 34 others, 63 trios. 215 SNPs.

gmr4M <- GetMaybeRel(GenoM = check4M,
                     Err = error_rate,
                     # SeqList = outfull,
                     Module = "par",
                     # MaxMismatch = NA,
                     Complex = "simp",
                     LifeHistData = LH_Data,
                     quiet = FALSE,
                     Tassign = 1.0,
                     # Tfilter = -100,
                     MaxPairs = 7*nrow(check4M))
# 135 PO pairs, 36 others, 60 trios. 182 SNPs.

gmr5M <- GetMaybeRel(GenoM = check5M,
                     Err = error_rate,
                     # SeqList = outfull,
                     Module = "par",
                     # MaxMismatch = NA,
                     Complex = "simp",
                     LifeHistData = LH_Data,
                     quiet = FALSE,
                     Tassign = 1.0,
                     # Tfilter = -100,
                     MaxPairs = 7*nrow(check5M))
# 124 PO pairs, 46 others, 64 trios. 154 SNPs.

##### and what happens when you increase the error? #####
# seems like it'll be pretty difficult to find sites where there is NO missing data, 
# especially as the pool of individuals increases (assigning ALL of the F1's). so 
# we're going to have to find a way to work with some missing data involved. let's
# first try to do this by slightly increasing the error to 0.05.

error_rate <- c(0.05, 0.05, 0.05)
gmr1M <- GetMaybeRel(GenoM = check1M,
                     Err = error_rate,
                     # SeqList = outfull,
                     Module = "par",
                     # MaxMismatch = NA,
                     Complex = "simp",
                     LifeHistData = LH_Data,
                     quiet = FALSE,
                     Tassign = 1.0,
                     # Tfilter = -100,
                     MaxPairs = 7*nrow(check1M))
# 150 PO pairs, 48 others, 90 trios. 437 SNPs.

gmr2M <- GetMaybeRel(GenoM = check2M,
                     Err = error_rate,
                     # SeqList = outfull,
                     Module = "par",
                     # MaxMismatch = NA,
                     Complex = "simp",
                     LifeHistData = LH_Data,
                     quiet = FALSE,
                     Tassign = 1.0,
                     # Tfilter = -100,
                     MaxPairs = 7*nrow(check2M))
# 149 PO pairs, 59 others, 93 trios. 279 SNPs.

gmr3M <- GetMaybeRel(GenoM = check3M,
                     Err = error_rate,
                     # SeqList = outfull,
                     Module = "par",
                     # MaxMismatch = NA,
                     Complex = "simp",
                     LifeHistData = LH_Data,
                     quiet = FALSE,
                     Tassign = 1.0,
                     # Tfilter = -100,
                     MaxPairs = 7*nrow(check3M))
# 137 PO pairs, 81 others, 93 trios. 215 SNPs.

gmr4M <- GetMaybeRel(GenoM = check4M,
                     Err = error_rate,
                     # SeqList = outfull,
                     Module = "par",
                     # MaxMismatch = NA,
                     Complex = "simp",
                     LifeHistData = LH_Data,
                     quiet = FALSE,
                     Tassign = 1.0,
                     # Tfilter = -100,
                     MaxPairs = 7*nrow(check4M))
# 136 PO pairs, 88 others, 100 trios. 182 SNPs.

gmr5M <- GetMaybeRel(GenoM = check5M,
                     Err = error_rate,
                     # SeqList = outfull,
                     Module = "par",
                     # MaxMismatch = NA,
                     Complex = "simp",
                     LifeHistData = LH_Data,
                     quiet = FALSE,
                     Tassign = 1.0,
                     # Tfilter = -100,
                     MaxPairs = 7*nrow(check5M))
# 131 PO pairs, 91 others, 101 trios. 154 SNPs.

# alright, here's the deal:
    # with no missing data and low error (1%), we didn't get to 95. (max 87)
    # with missing data and low error (1%), we didn't get to 95 (max 65)
    # with missing data and high error (5%), we got to 101! (max 101)

# now the question becomes:
    # 1. with no missing data, how high does the error have to be to get enough trios?
    # 2. can we use a more strict min depth filter to keep higher quality loci?
    # 3. will more loci be able to assign the 95 trios with low specified error?
    # 4. this MAY require reducing the thinning. what is the MINIMUM thinning I should use?

# 1. with no missing data, how high does the error have to be to get enough trios?
# with no missing data and an error of 0.01, we got the first line's assignments.
# second line's assignments are with 0.05 error.
# let's split the difference and try with 0.025 error. (third line)
error_rate <- c(0.025, 0.025, 0.025)
gmr1M_nomiss <- GetMaybeRel(GenoM = check1M_nomiss,
                            Err = error_rate,
                            # SeqList = outfull,
                            Module = "par",
                            # MaxMismatch = NA,
                            Complex = "simp",
                            LifeHistData = LH_Data,
                            quiet = FALSE,
                            Tassign = 1.0,
                            # Tfilter = -100,
                            MaxPairs = 7*nrow(check1M_nomiss))
# 157 PO pairs, 34 others, 80 trios. 129 SNPs.
# 161 PO pairs, 99 others, 112 trios. 129 SNPs.
# 153 PO pairs, 52 others, 86 trios. 129 SNPs.

gmr2M_nomiss <- GetMaybeRel(GenoM = check2M_nomiss,
                            Err = error_rate,
                            # SeqList = outfull,
                            Module = "par",
                            # MaxMismatch = NA,
                            Complex = "simp",
                            LifeHistData = LH_Data,
                            quiet = FALSE,
                            Tassign = 1.0,
                            # Tfilter = -100,
                            MaxPairs = 7*nrow(check2M_nomiss))
# 166 PO pairs, 58 others, 87 trios. 79 SNPs
# 177 PO pairs, 132 others, 143 trios. 79 SNPs
# 180 PO pairs, 81 others, 106 trios. 79 SNPs.

gmr3M_nomiss <- GetMaybeRel(GenoM = check3M_nomiss,
                            Err = error_rate,
                            # SeqList = outfull,
                            Module = "par",
                            # MaxMismatch = NA,
                            Complex = "simp",
                            LifeHistData = LH_Data,
                            quiet = FALSE,
                            Tassign = 1.0,
                            # Tfilter = -100,
                            MaxPairs = 7*nrow(check3M_nomiss))
# 166 PO pairs, 77 others, 85 trios. 65 SNPs
# 168 PO pairs, 129 others, 129 trios. 65 SNPs.
# 175 PO pairs, 103 others, 104 trios. 65 SNPs.

gmr4M_nomiss <- GetMaybeRel(GenoM = check4M_nomiss,
                            Err = error_rate,
                            # SeqList = outfull,
                            Module = "par",
                            # MaxMismatch = NA,
                            Complex = "simp",
                            LifeHistData = LH_Data,
                            quiet = FALSE,
                            Tassign = 1.0,
                            # Tfilter = -100,
                            MaxPairs = 7*nrow(check4M_nomiss))
# 166 PO pairs, 92 others, 80 trios. 51 SNPs.
# 137 PO pairs, 140 others, 124 trios. 51 SNPs.
# 165 PO pairs, 118 others, 109 trios. 51 SNPs.

gmr5M_nomiss <- GetMaybeRel(GenoM = check5M_nomiss,
                            Err = error_rate,
                            # SeqList = outfull,
                            Module = "par",
                            # MaxMismatch = NA,
                            Complex = "simp",
                            LifeHistData = LH_Data,
                            quiet = FALSE,
                            Tassign = 1.0,
                            # Tfilter = -100,
                            MaxPairs = 7*nrow(check5M_nomiss))
# 170 PO pairs, 100 others, 63 trios. 43 SNPs.
# 113 PO pairs, 128 others, 89 trios. 43 SNPs.
# 162 PO pairs, 112 others, 79 trios. 43 SNPs.








