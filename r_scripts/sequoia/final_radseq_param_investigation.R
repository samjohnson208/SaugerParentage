# this code is from Sequoia_ContamFilt_mindep8_md5_RADseqErr_ALLINDS.R
# but i need to pull in the new crosses that we discovered and updated
# on 0210 and 021126. 

# the new cross lookup table is temporarily stored here:
#/Users/samjohnson/Desktop/sar_2015_filt_split_pair.csv

# UPDATE WITH PERMANENT LOCATION: /Users/samjohnson/Documents/GeneticData/F0_CROSSES_021326


# we'll start with the same genotype matrix and LH data that we've been working
# with so far, but we'll have to filter to only the correct F0's in the new lookup.

# save.image(file = "resultsplotted_EOD_020426.RData")
# this ^ is from extrapolation_120925.R, the script where we figured out the two
# pronged sequoia analysis and where plotting code is for the results of this proj.

# but that only has the thin100K matrix, so we need all these

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

# these are read in from Sequoia_ContamFilt_mindep8.R, object names changed from including min8 to not.

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
dim(LH_Test) #(comes from resultsplotted_EOD_020426.RData)

# filter the genotype matrices so that it only includes the ids from LH_Test$ID
mat_thin25K <- mat_thin25K[rownames(mat_thin25K) %in% LH_Test$ID, , drop = FALSE]
mat_thin50K <- mat_thin50K[rownames(mat_thin50K) %in% LH_Test$ID, , drop = FALSE]
mat_thin75K <- mat_thin75K[rownames(mat_thin75K) %in% LH_Test$ID, , drop = FALSE]
mat_thin100K <- mat_thin100K[rownames(mat_thin100K) %in% LH_Test$ID, , drop = FALSE]
mat_thin200K <- mat_thin200K[rownames(mat_thin200K) %in% LH_Test$ID, , drop = FALSE]
mat_thin300K <- mat_thin300K[rownames(mat_thin300K) %in% LH_Test$ID, , drop = FALSE]
mat_thin400K <- mat_thin400K[rownames(mat_thin400K) %in% LH_Test$ID, , drop = FALSE]
mat_thin500K <- mat_thin500K[rownames(mat_thin500K) %in% LH_Test$ID, , drop = FALSE]
mat_thin600K <- mat_thin600K[rownames(mat_thin600K) %in% LH_Test$ID, , drop = FALSE]
mat_thin700K <- mat_thin700K[rownames(mat_thin700K) %in% LH_Test$ID, , drop = FALSE]
mat_thin800K <- mat_thin800K[rownames(mat_thin800K) %in% LH_Test$ID, , drop = FALSE]
mat_thin900K <- mat_thin900K[rownames(mat_thin900K) %in% LH_Test$ID, , drop = FALSE]
mat_thin1M <- mat_thin1M[rownames(mat_thin1M) %in% LH_Test$ID, , drop = FALSE]
mat_thin2M <- mat_thin2M[rownames(mat_thin2M) %in% LH_Test$ID, , drop = FALSE]
mat_thin3M <- mat_thin3M[rownames(mat_thin3M) %in% LH_Test$ID, , drop = FALSE]
mat_thin4M <- mat_thin4M[rownames(mat_thin4M) %in% LH_Test$ID, , drop = FALSE]
mat_thin5M <- mat_thin5M[rownames(mat_thin5M) %in% LH_Test$ID, , drop = FALSE]

# 206 indivs # 112 F0 2015 parents spawned, 94 test F1 offspring
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

##### ---- HWE Filtering and Final GenoMat Maintenence ---- #####
# except no HWE filtering here.

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

dim(mat_thin25K) # 1100 <- 1097
dim(mat_thin50K) # 1041 <- 1038
dim(mat_thin75K) # 1005 <- 1002
dim(mat_thin100K) # 943 <- 940
dim(mat_thin200K) # 808 <- 805
dim(mat_thin300K) # 711 <- 708
dim(mat_thin400K) # 648 <- 646
dim(mat_thin500K) # 594 <- 592
dim(mat_thin600K) # 546 <- 543
dim(mat_thin700K) # 509 <- 506
dim(mat_thin800K) # 476 <- 475
dim(mat_thin900K) # 446 <- 444
dim(mat_thin1M) # 418 <- 417
dim(mat_thin2M) # 273 <- 272
dim(mat_thin3M) # 211 <- 210
dim(mat_thin4M) # 175 <- 175
dim(mat_thin5M) # 151 <- 151


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

dim(check_thin25K)
dim(check_thin50K)
dim(check_thin75K)
dim(check_thin100K)
dim(check_thin200K)
dim(check_thin300K)
dim(check_thin400K)
dim(check_thin500K)
dim(check_thin600K) 
dim(check_thin700K) 
dim(check_thin800K) 
dim(check_thin900K) 
dim(check_thin1M) 
dim(check_thin2M) 
dim(check_thin3M) 
dim(check_thin4M)
dim(check_thin5M)

##### ---- RADseq Error Matrix Exploration ---- #####
setwd("/Users/samjohnson/Desktop/")
cross_lookup_new <- read.csv(file = "sar_2015_filt_split_pair.csv", header = TRUE)
cross_lookup_new <- cross_lookup_new %>% 
  select(Sample_ID, Crossed_with, Pair) %>% 
  rename(pair = Pair)

library(sequoia)

# need this just to establish the ageprior, don't worry about error here.
seq_test_params <- sequoia(GenoM = check_thin100K,
                           LifeHistData = LH_Test,
                           Module = "ped",
                           Complex = "full",
                           Herm = "no",
                           # Err = errM,
                           UseAge = "yes",
                           args.AP=list(Discrete = TRUE, 
                                        MinAgeParent = 1, MaxAgeParent = 1),
                           CalcLLR = TRUE,
                           StrictGenoCheck = TRUE,
                           DummyPrefix = c("F", "M"),
                           Tfilter = -2,
                           Tassign = 0.5)

# create vectors of e0 and e1 values
e0_vals <- c(0.005, 0.01, 0.025, 0.05, 0.075, 0.1)
e1_vals <- c(0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5)

# create empty dfs to display
assign_rate_mat <- matrix(NA, nrow = length(e1_vals), ncol = length(e0_vals),
                          dimnames = list(paste0("e1_", e1_vals),
                                          paste0("e0_", e0_vals)))
acc_rate_mat <- assign_rate_mat  # same structure

# total number of offspring
n_offspring <- 95

# read in cross lookup table
cross_lookup_new <- read.csv(file = "sar_2015_filt_split_pair.csv", header = TRUE)
cross_lookup_new <- cross_lookup_new %>% 
  select(Sample_ID, Crossed_with, Pair) %>% 
  rename(pair = Pair)

# nested loop for all combinations:
for (i in seq_along(e1_vals)) {
  for (j in seq_along(e0_vals)) {
    
    e1_val <- e1_vals[i]
    e0_val <- e0_vals[j]
    
    # print text to display each combination, track progress
    cat("Running GetMaybeRel for E0 =", e0_val, "and E1 =", e1_val, "\n")
    
    # create the error matrix for that combination
    errM <- Err_RADseq(E0 = e0_val, E1 = e1_val, Return = 'matrix')
    
    # run gmr (change dataset here only! how great is that?!)
    gmr <- GetMaybeRel(GenoM = check_thin5M,
                       Err = errM,
                       Module = "par",
                       Complex = "simp",
                       LifeHistData = LH_Test,
                       AgePrior = seq_test_params[["AgePriors"]],
                       quiet = TRUE,
                       Tassign = 0.5,
                       Tfilter = -2,
                       Herm = "no",
                       MaxPairs = 7 * nrow(check_thin5M))
    
    # if there are no trios found, place 0 in that place for assign, NA for accuracy
    if (is.null(gmr[["MaybeTrio"]]) || nrow(gmr[["MaybeTrio"]]) == 0) {
      cat("   No trios found for this combination. Skipping...\n")
      assign_rate_mat[i, j] <- 0
      acc_rate_mat[i, j] <- NA
      next
    }
    
    # how many unique test f1's were assigned to parents?
    assigned_ids <- unique(gmr[["MaybeTrio"]]$id[grepl("^SAR_15_67", gmr[["MaybeTrio"]]$id)])
    n_assigned <- length(assigned_ids) # define it as an object
    
    # use that object to calculate assignment rate
    assign_rate <- (n_assigned / n_offspring) * 100 
    
    # make trios an object, create the pair entry, check it against the lookup...
    trios <- gmr[["MaybeTrio"]] %>%
      mutate(pair = paste(pmin(parent1, parent2), pmax(parent1, parent2), sep = "_")) %>%
      left_join(cross_lookup %>%
                  select(pair) %>%
                  distinct() %>%
                  mutate(valid_cross = TRUE),
                by = "pair") %>%
      mutate(valid_cross = ifelse(is.na(valid_cross), FALSE, valid_cross)) %>%
      # ...and filter out any placeholders. 
      filter(!LLRparent1 %in% c(555, -555),
             !LLRparent2 %in% c(555, -555),
             !LLRpair    %in% c(555, -555))
    
    # if we have assigned individuals, grab the number of valid crosses and calculate
    # the accuracy rate by dividing it by the number of unique test f1's that got assigned
    if (n_assigned > 0) {
      valid_ids <- unique(trios$id[trios$valid_cross])
      n_valid <- length(valid_ids)
      acc_rate <- (n_valid / n_assigned) * 100
    } else {
      acc_rate <- NA
    }
    
    # place those values in the correct cell of each summary table.
    assign_rate_mat[i, j] <- assign_rate
    acc_rate_mat[i, j] <- acc_rate
  }
}


# print results, round to three decimal places
cat("\nAssignment rate matrix (%):\n")
print(round(assign_rate_mat, 3))

cat("\nAccuracy rate matrix (%):\n")
print(round(acc_rate_mat, 3))

# create composite score matrix
composite_mat <- (assign_rate_mat * acc_rate_mat) / 100

# use the same row and column names as the other two result tables
dimnames(composite_mat) <- dimnames(assign_rate_mat)

# print results, round to three decimal places
cat("\nComposite score matrix (%):\n")
print(round(composite_mat, 3))

# notes:


##### ---- ---- #####






