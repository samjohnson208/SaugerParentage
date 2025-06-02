### --- hiphop_on_simdata.R --- ###
### created by SPJ on 060225 ###
# purpose: check out hiphop on simulated data for f0-f1 and f0-f2

setwd("/Users/samjohnson/Documents/Sauger_042225/SaugerParentage/sims_hiphop")
genotypes <- read.csv(file = "all_inds_matrix.csv", header = TRUE)
rownames(genotypes) <- genotypes[,1]
genotypes <- genotypes[, 2:ncol(genotypes)]
