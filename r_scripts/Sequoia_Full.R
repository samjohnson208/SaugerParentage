# unsure if this will work first try. worth a shot.

# install.packages("sequoia")
# 65
library(sequoia)

# install.packages("dplyr")
# 65
library(dplyr)

setwd("/project/ysctrout/hatchsauger/sam_sai_svit/Sequoia_Inp/maf30_miss9")
#setwd("/Users/samjohnson/Desktop/")

# read in genotype matrix
mat <- read.table(file = "variants_maf30_miss9.012_conv", header = FALSE, sep = "\t", 
                  na.strings = c("NA", "-1"))

# correct row names as sample id's
mat<- mat[, -1]
gmmat <- as.matrix(mat)
ind <- read.table(file = "variants_maf30_miss9.012.indv", header = FALSE)
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

# generate snp stats and plots for the checked genotype matrix and output to pdf
pdf("plots.pdf")
stats <- SnpStats(gmmat_po_check, Pedigree = NULL, Duplicates = NULL, Plot = TRUE, quiet = FALSE)
dev.off()

# Not conditioning on any pedigree^

# run sequoia on the checked genotype matrix
outfull <- sequoia(GenoM = gmmat_po_check, Module = 'ped', StrictGenoCheck = TRUE, CalcLLR = TRUE)

# run GetMaybeRel() on the sequoia output
# gmr <- GetMaybeRel(outfull, GenoM = gmmat)
gmr <- GetMaybeRel(GenoM = gmmat_po_check)

# save output to current wd 
save(outfull, file = "Sequoia_OutFull_032425.RData")
save(gmr, file = "Sequoia_GetMayRel_032425.RData")

### Exploring Output ### 
# setwd("/Users/samjohnson/Desktop")
# load("Sequoia_GetMayRel_022725.RData")
# load("Sequoia_OutFull_022725.RData")
# 
# library(sequoia)
# SummarySeq(outfull)
# SummarySeq(outfull)

# See example code in ?GetMaybeRel()

# setwd("/Users/samjohnson/Desktop/")
# load("Sequoia_GetMayRel_031325.RData")
# load("Sequoia_OutFull_031325.RData")


setwd("/Users/samjohnson/Desktop/")
load("Sequoia_GetMayRel_032425.RData")
load("Sequoia_OutFull_032425.RData")


