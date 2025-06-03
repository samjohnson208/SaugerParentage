### --- hiphop_on_simdata.R --- ###
### created by SPJ on 060225 ###
# purpose: check out hiphop on simulated data for f0-f1 and f0-f2

setwd("/Users/samjohnson/Documents/Sauger_042225/SaugerParentage/sims_hiphop")
genotypes <- read.csv(file = "all_inds_matrix.csv", header = TRUE)
rownames(genotypes) <- genotypes[,1]
genotypes <- genotypes[, 2:ncol(genotypes)]

true_parents <- read.table(file = "true_parents_1000_1.txt", header = FALSE)

# where i stand right now
# have the genotypes in, as well as a record of true parents for both generations
# for 1000 loci, and all parents sampled

# what this exercise will do (next steps)
# the first thing to do is to make the individual tables for each run:
# f1 to f0 (hatch and wild, remember that cross prop is 0.2)
# f2 to f0

f1f0ind <- read.csv(file = "F1F0indivs.csv", header = TRUE)
f2f0ind <- read.csv(file = "F2F0indivs.csv", header = TRUE)


# then filter the genotype matrix to include indivs from those indiv tables

f1f0geno <- genotypes %>% 
    filter(rownames(genotypes) %in% f1f0ind$individual)
f2f0geno <- genotypes %>% 
  filter(rownames(genotypes) %in% f2f0ind$individual)

# then run hothiphop on both generations individually (f1 to f0
# and f2 to f0) all going to be false.

# then we'll need to plot all of those scores and color by valid cross
# five categories (violins) to plot are going to be (f1 to f0, valid cross = T/F,
# and f2 to f0 (all F), f2 to f1 valid cross = T/F)

# what this exercise does NOT include
# situations where f0 x f0 made f2 individuals, or cases where f0 x f1 made 
# f2 individuals. we will likely need to understand what those scores look like.

# then after ALL OF THIS ^ 
# we'll need to analyze the effects of PTPS on these relationships and we'll start
# pulling some of the parents in differing frequencies and replicating.


# final thoughts 6/2
# hey look. don't try to do too much. our questions are this:
# 1. Does hiphop work on the simulated data? (run it and valid cross)
# 2. What does our resolution look like when we have various ptps? (shoot josh a message tomorrow morning?)
# 3. What happens when we try to assign F0's to f

# DON'T FORGET TO PUSH

