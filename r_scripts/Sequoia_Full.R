# unsure if this will work first try. worth a shot.

# install.packages("sequoia")
# 65
library(sequoia)

# install.packages("dplyr")
# 65
library(dplyr)

setwd("/project/ysctrout/hatchsauger/sam_sai_svit/Sequoia_Inp")

# read in genotype matrix
mat <- read.table(file = "variants_maf1_miss9.012_conv", header = FALSE, sep = "\t", 
                  na.strings = c("NA", "-1"))

# correct row names as sample id's
mat<- mat[, -1]
gmmat <- as.matrix(mat)
ind <- read.table(file = "variants_maf1_miss9.012.indv", header = FALSE)
str(ind)
ind <- ind %>% 
  rename(sample = V1)
rownames(gmmat) <- ind$sample

# check genotype matrix for samples/loci to be excluded
check <- CheckGeno(gmmat, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))    
# There are 2 individuals scored for <5% of SNPs, these WILL BE IGNORED
# In addition, there are 6 individuals scored for <20% of SNPs,it is advised to treat their assignments with caution
# After exclusion, There are  1182  individuals and  6468  SNPs.

# generate snp stats and plots for the genotype matrix and output to pdf
pdf("plots.pdf")
stats <- SnpStats(gmmat, Pedigree = NULL, Duplicates = NULL, Plot = TRUE, quiet = FALSE)
dev.off()

# Not conditioning on any pedigree

# run sequoia on the genotype matrix
outfull <- sequoia(GenoM = gmmat, Module = 'ped', MaxSibIter = 42, StrictGenoCheck = TRUE, CalcLLR = TRUE)

# run GetMaybeRel() on the sequoia output
# gmr <- GetMaybeRel(outfull, GenoM = gmmat)
gmr <- GetMaybeRel(GenoM = gmmat)

# save output to current wd 
save(outfull, file = "Sequoia_OutFull_022725.RData")
save(gmr, file = "Sequoia_GetMayRel_022725.RData")

### Exploring Output ### 
setwd("/Users/samjohnson/Desktop")
load("Sequoia_GetMayRel_022725.RData")
load("Sequoia_OutFull_022725.RData")

library(sequoia)
SummarySeq(outfull)
SummarySeq(outfull)

SnpStats()

# hm... don't get it. will have to come back.