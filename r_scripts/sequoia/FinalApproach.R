##### ---- FinalApproach.R by SPJ 022426 #####
## PURPOSE: over the past few weeks, i have come up with this two-pronged approach 
# that relies on 1. CalcPairLL() -> LLtoProb() and 2. sequoia() -> GetMaybeRel()
# -> GetRelM() to 1. come up with probabilities that each pair of inds is related
# according to each possible relationship type, and 2. infer the most likely pedigree
# and generate a relatedness matrix from that. (see extrapolation_120925.R)

# because it'll be quite difficult to articulate that approach and to interpret
# results in a talk or paper format, we've decide it's best to combine the two approaches,
# or to pick one of the two. it appears as though you can run CalcPairLL() conditional
# on an inferred pedigree, so the plan is to do that on each of the groups and then
# make the same plots as before. number of PO assignments, 

##### ---- packages #####
library(sequoia)
library(tidyverse)
library(ggplot2)
library(ggpattern)

##### ---- load and tailor genotype matrix #####
mat_thin100K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012_conv",
                           header = FALSE, sep = "\t", na.strings = c("NA", "-1"))

mat_thin100K <- mat_thin100K[, -1]

mat_thin100K <- as.matrix(mat_thin100K)

dim(mat_thin100K)

ind100K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012.indv", header = FALSE)

ind100K <-  ind100K %>% 
  rename(sample = V1)

rownames(mat_thin100K) <- ind100K$sample

pos100K <- read.table(file = "rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.012.pos", header = FALSE)

pos100K$position <-  paste(pos100K$V1, pos100K$V2, sep = "_")

colnames(mat_thin100K) <- pos100K$position







