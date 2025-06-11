## Sequoia_ReducedNLoci.R by SPJ 060225
## PURPOSE: to understand problems assoc. with Sequoia. first step is to attempt
          # running sequoia with reduced numbers of loci. i am concerned that the
          # problems i have been having are due to the use of too many loci given
          # what sequoia specifies in its documentation.
## USAGE: Rscript Sequoia_ReducedNLoci.R

# install.packages("sequoia")
# 65
library(sequoia)

# install.packages("dplyr")
# 65
library(dplyr)

# going to run it locally first. then cluster if run time is too long.
setwd("/Users/samjohnson/Desktop/")

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
outfull <- sequoia(GenoM = gmmat_po_check, Module = 'ped', StrictGenoCheck = TRUE, CalcLLR = TRUE)
gmr <- GetMaybeRel(GenoM = gmmat_po_check)

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
  MaxMismatch = 100,   # <- RELAX threshold here
  quiet = FALSE,
  Tassign = 0, 
  Tfilter = -10,
  MaxPairs = 7*nrow(filtered_gmmat_po_check)
)
# got only 46 po pairs that time...


#let's try 10%
samp_loci <- sample(ncol(filtered_gmmat_po_check), 129)
filtered_gmmat_po_check_10 <- filtered_gmmat_po_check[, samp_loci]
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
