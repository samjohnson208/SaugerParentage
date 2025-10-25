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

table(LH_All$BirthYear)
# 263 F0's (15 and 16), 6436 was not spawned, so it's not included here.
# 334 F1's (juveniles from fall 15 and spawn agg fish from spring 21)
# NO TEST F1's included
# 480 F2's (19, 20, 21, all over)
# = 1077
# note: will's samples excluded, possible hybrids excluded, WF fish excluded
# THIS IS READY NOW FOR GETMAYBEREL

# filter GenoM so it only includes ids from LH_All$ID
dim(mat_thin100K) # 1184, 943
mat_thin100K <- mat_thin100K[rownames(mat_thin100K) %in% LH_All$ID, , drop = FALSE]
dim(mat_thin100K) # 1060 943
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
    filter(ID %in% rownames(mat_thin100K))
dim(mat_thin100K) # 1060 943
dim(LH_All)       # 1090 3
# okay, that's remedied. great.


##### ----- ---- #####

##### ---- HWE Filtering and Final GenoMat Maintenence ---- #####
# except no HWE filtering here.

# check for all heterozygous sites since those likely will not be informative.
all_ones_100 <- apply(mat_thin100K, 2, function(col) all(col == 1))
# remove all het. sites
mat_thin100K <- mat_thin100K[, !all_ones_100]

dim(mat_thin100K) # retained all

check_thin100K <- CheckGeno(mat_thin100K, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))
# ✖ There are 2 individuals scored for <5% of SNPs, these WILL BE IGNORED
# ✖ In addition, there are 6 individuals scored for <20% of SNPs, it is advised to treat their assignments with caution
# ℹ After exclusion, There are  1058  individuals and  943  SNPs.

# should consider filtering for per individual md!!!
# let's just let it fly right now and see what happens. we'll tighten that up later.


##### ----- ---- #####

##### ---- Missing Data per Individual ---- #####
# first runs say that there are six individuals genotyped for <20% of snps. 
# let's remedy that.

miss_per_ind <- read.table(file = "miss_per_indv_rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.imiss", header = TRUE)

inds_to_keep <- miss_per_ind$INDV[miss_per_ind$F_MISS < 0.2]
inds_to_keep # individuals who are genotyped for 80% or samples or MORE
length(inds_to_keep) # 1154
dim(miss_per_ind) # 1184, means 30 inds are filtered out...

dim(check_thin100K) # 1058 samples from JUST the groups we're interested in (F0, F1Spawn, F1Juv, F2)
check_thin100K <- check_thin100K[rownames(check_thin100K) %in% inds_to_keep, ]
dim(check_thin100K) # 1030 samples remaining.

##### ---- ---- #####

##### ---- sequoia() and GetMaybeRel() ---- #####

errM <- Err_RADseq(E0 = 0.075, E1 = 0.025, Return = 'matrix')
errM

seq <- sequoia(GenoM = check_thin100K, 
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
# run two after removing individuals with 20% or more missing data. surprised that
# that there were fewer assignments, though there are fewer inds i suppose...
# I hope that this also means that more of these assignments are accurate.


# overlapping gens? going to have to fix this. bet that ups the number of assignments.
# that's obviously the next step.

gmr <- GetMaybeRel(GenoM = check_thin100K,
                   Err = errM,
                   Module = "ped",
                   Complex = "full",
                   LifeHistData = LH_All,
                   quiet = TRUE,
                   Tfilter = -2,
                   Tassign = 1.0,
                   MaxPairs = 7 * nrow(check_thin1M))


##### ----- ---- #####

# notes:
# work on generation times and making the ageprior more informative
# saying you get FS at agedif = 2. ask Jisca? just try some more stuff out first.


# could try subsetting to just F0 and F1, seeing if we can extract more relationships,
# similar to what we were doing with the test set? do pairwise generations?


# duplicate check? see the seq objects (nothing to sweat. just some inds that are similar)
# would expect them to be siblings or something. this is a low div pop.

# per-sample md? DONE
# 6 indivs scored for < 20 % of snps
# 2 indivs genotyped for < 5 % of snps
# should be fixed now.


