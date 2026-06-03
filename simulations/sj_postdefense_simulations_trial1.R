##### ---- sj_postdefense_simulations_trial1.R ---- #####
### by SPJ 060226 ###

library(tidyverse)

##### ---- notes on goals and approach ---- #####

# i have used code from JPJ in SaugerParentage/r_scripts/sfs/sfs_all.R to estimate
# the alpha and beta parameters of the site frequency spectrum for Will's above-dam
# sauger. here, i will use those alpha and beta params to create a distrubition, and
# draw alleles from that distribution using the binomial distribution. i will create
# an F0 generation, cross everyone in that F0 generation to create an F1 group.
# then i'll take 95 of those F1s and all of their parents, and that's the test group.
# see how sequoia does (best case), then we start pulling some parents (dropping ptps),
# and adding some false ones (100% ptps + some false).

# for the real F0-F1 PO scenario, you have some unknown ptps, but no extras. that's
# the thing. so i'm not sure how applicable the extras are. they are for the F1-F2
# scenario.

# true scenarios for the data...
# F0-F1: 100% ptps for the stocked fish, for the wild fish, you have some unknown ptps.
    # so for the whole F1 group, you have some unknown ptps.
# F1-F2: you have BOTH less than 100% ptps, AND some extra false parents

### scenarios we're testing ###
# 1. best case scenario, 100% ptps with no extras (test)
# 2. some ptps < 100% with no extras (replicating PO for F0-F1)
# 3. 100% ptps plus some extras (THIS IS THE TEST GROUP, BECAUSE THE TEST F1s
   # came from only a portion of all the crosses that were made. 60/110 parents
   # were inferred as parents of the test offspring.)
# 4. also need less than ideal ptps WITH some extras (for F0-F1 PO)

# THEN WE NEED TO DO THIS SAME THING WITH THE GPs AND SEE HOW IT PERFORMS.

##### ----  ---- #####

set.seed(051999)

##### ---- params for sfs to simulate F0 generation ---- #####

# these were drawn from the empirical sfs for will's above-dam fish. i took the 
# frequency of the alternate allele, and fit the beta dist to that. wondering if
# that's the right approach... not quite sure.

alpha <- 0.0657024
beta <- 0.2714633

##### ----  ---- #####

##### ---- simulate F0 generation ---- #####

# specify params here
nloci <- 10000
f0_inds <- 500
f1_inds <- 500

# establish matrix of true parent id's
true_parents <- data.frame(
  offspring = rep(NA_character_, f0_inds + f1_inds),
  parent1   = rep(NA_character_, f0_inds + f1_inds),
  parent2   = rep(NA_character_, f0_inds + f1_inds),
  stringsAsFactors = FALSE
)


# establish genotype matrix
geno_mat_f0 <- matrix(NA, f0_inds, nloci)
for (i in 1:nloci) {
  freq <- rbeta(1, alpha, beta)  # beta estimates come from wcr data
  for (j in 1:f0_inds) {
    geno_mat_f0[j, i] <- sum(rbinom(2, 1, freq))
  }
}

true_parents$offspring[1:f0_inds] <- paste0("f0_", 1:f0_inds)
rownames(geno_mat_f0) <- true_parents$offspring[1:f0_inds]

# specify sexes of the parents
sex <- ifelse(1:f0_inds %% 2 == 0, "M", "F")
true_parents$sex <- NA
true_parents$sex[1:f0_inds] <- sex

f0_males <- true_parents %>% 
    filter(sex == "M")
f0_males <- f0_males$offspring

f0_females <- true_parents %>% 
    filter(sex == "F")
f0_females <- f0_females$offspring

##### ---- ---- #####

##### ---- cross F0 generation (create F1 generation) ---- #####

# assign sample names to f1 individuals
true_parents$offspring[(f0_inds+1):(f0_inds+f1_inds)] <- paste0("f1_", 1:f1_inds)

# sample F0 parents and place their names into the true_parents table 
for(i in 1:f1_inds){
  
  father <- sample(f0_males, 1)
  mother <- sample(f0_females, 1)
  
  row <- f0_inds + i
  
  true_parents$parent1[row] <- father
  true_parents$parent2[row] <- mother
}

# let's make sure that we have a cross column that describes the true parents
# and sorts them numerically.
true_parents$parent_cross <- NA

for(i in (f0_inds+1):(f0_inds+f1_inds)){
  
  p1 <- true_parents$parent1[i]
  p2 <- true_parents$parent2[i]
  
  id1 <- as.numeric(sub("f0_", "", p1))
  id2 <- as.numeric(sub("f0_", "", p2))
  
  true_parents$parent_cross[i] <- paste(
    paste0("f0_", pmin(id1, id2)),
    paste0("f0_", pmax(id1, id2)),
    sep="__"
  )
}

id_to_row <- setNames(1:f0_inds, paste0("f0_", 1:f0_inds))

### SIMULATE F1 GENOTYPES FROM F0 PARENTS ###
geno_mat_f1 <- matrix(NA, f1_inds, nloci)

rownames(geno_mat_f1) <- true_parents$offspring[(f0_inds + 1):(f0_inds + f1_inds)]

for (i in 1:f1_inds){
  # grab the parents from true_parents
  row <- f0_inds + i
  
  # select the parents from true_parents
  p1 <- true_parents[row, "parent1"]
  p2 <- true_parents[row, "parent2"]
  
  # make the crosses (shoutout jpj's mendel skills)
    for(j in 1:nloci){
      
      g1 <- geno_mat_f0[id_to_row[p1], j]
      g2 <- geno_mat_f0[id_to_row[p2], j]
      
      if      (geno_mat_f0[p1,j]==0 && geno_mat_f0[p2,j]==0) { geno_mat_f1[i,j] <- 0}
      else if (geno_mat_f0[p1,j]==0 && geno_mat_f0[p2,j]==1) { geno_mat_f1[i,j] <- rbinom(1,1,0.5) }
      else if (geno_mat_f0[p1,j]==0 && geno_mat_f0[p2,j]==2) { geno_mat_f1[i,j] <- 1 }
      else if (geno_mat_f0[p1,j]==1 && geno_mat_f0[p2,j]==0) { geno_mat_f1[i,j] <- rbinom(1,1,0.5) }
      else if (geno_mat_f0[p1,j]==1 && geno_mat_f0[p2,j]==1) { geno_mat_f1[i,j] <- sum(rbinom(2,1,0.5)) }
      else if (geno_mat_f0[p1,j]==1 && geno_mat_f0[p2,j]==2) { geno_mat_f1[i,j] <- rbinom(1,1,0.5) + 1 }
      else if (geno_mat_f0[p1,j]==2 && geno_mat_f0[p2,j]==0) { geno_mat_f1[i,j] <- 1 }
      else if (geno_mat_f0[p1,j]==2 && geno_mat_f0[p2,j]==1) { geno_mat_f1[i,j] <- rbinom(1,1,0.5) + 1 }
      else if (geno_mat_f0[p1,j]==2 && geno_mat_f0[p2,j]==2) { geno_mat_f1[i,j] <- 2 }
      }
  }

geno_mat_f0f1 <- rbind(geno_mat_f0, geno_mat_f1)
dim(geno_mat_f0f1)

##### ---- ---- #####

dim(geno_mat_f0f1)

table(table(unique(true_parents$parent1)))
# alright, so we've gotten to each of the generations, but sampling a group of 
# 95 f1's and all of their parents would not be at all reminiscent of the test
# group, because we'd have almost all 500 f0 individuals in that matrix. why don't
# we take a subset of the F0 population, cross them together by randomly drawing
# pairs, and i imagine at that point we'll have some duplicates. then we can make
# test offspring, some of those will be full or half sibs (good), and then we'll 
# be able to say which f0's were involved in the crosses, and which weren't. then
# we'll be pulling true ones, and we'll actually be able to add the false ones 
# back in since we'll have already created them...

##### ---- subsample the F0 generation to create broodstock families ---- #####
# sample 53 females, 59 males from the f0 generation 
# (replace = false to get unique inds)
test_f0males <- sample(f0_males, 59, replace = FALSE)
test_f0females <- sample(f0_females, 53, replace = FALSE)
test_f0parents <- c(test_f0males, test_f0females)

# create genotype matrix of broodstock contributors
geno_mat_testparents <- geno_mat_f0[rownames(geno_mat_f0) %in% test_f0parents, ]
dim(geno_mat_testparents)

# create the 61 unique crosses
test_crosses <- data.frame(female = sample(test_f0females, 61, replace = TRUE),
                           male = sample(test_f0males, 61, replace = TRUE),
                           stringsAsFactors = FALSE)

# how often were the parents used?
table(table(test_crosses$female))
table(table(test_crosses$male))
# okay... not entirely accurate to our situation, but let's push on.

##### ---- ---- #####

##### ---- create the test offspring ---- #####
# which of the 61 possible families will the 95 test offspring be from?
family_draws <- sample(1:61, 95, replace = TRUE)

# now go into the true_parents (pedigree) dataframe, add the 95 test offspring ids
true_parents$offspring[1001:1095] <- paste0("f1test_", 1:95)

# then for each one...
for(i in 1:95){
  row <- 1000 + i # define the individual that we're focusing on
  
  family <- family_draws[i] # grab that entry from the families vector (95 long)
  
  # then store the parents in the true_parents dataframe
  true_parents$parent1[row] <- test_crosses$male[family]
  true_parents$parent2[row] <- test_crosses$female[family]
  
  p1 <- true_parents$parent1[row]
  p2 <- true_parents$parent2[row]
  
  id1 <- as.numeric(sub("f0_", "", p1))
  id2 <- as.numeric(sub("f0_", "", p2))
  
  true_parents$parent_cross[row] <- paste(
    paste0("f0_", pmin(id1, id2)),
    paste0("f0_", pmax(id1, id2)),
    sep="__")
}

geno_mat_f1test <- matrix(NA, 95, nloci)
rownames(geno_mat_f1test) <- true_parents$offspring[1001:1095]

### simulate genotypes for the test f1 group
for(i in 1:95){
  row <- 1000 + i
  
  p1 <- true_parents$parent1[row]
  p2 <- true_parents$parent2[row]
  
  for(j in 1:nloci){
    
    if      (geno_mat_f0[p1,j]==0 && geno_mat_f0[p2,j]==0) { geno_mat_f1test[i,j] <- 0}
    else if (geno_mat_f0[p1,j]==0 && geno_mat_f0[p2,j]==1) { geno_mat_f1test[i,j] <- rbinom(1,1,0.5) }
    else if (geno_mat_f0[p1,j]==0 && geno_mat_f0[p2,j]==2) { geno_mat_f1test[i,j] <- 1 }
    else if (geno_mat_f0[p1,j]==1 && geno_mat_f0[p2,j]==0) { geno_mat_f1test[i,j] <- rbinom(1,1,0.5) }
    else if (geno_mat_f0[p1,j]==1 && geno_mat_f0[p2,j]==1) { geno_mat_f1test[i,j] <- sum(rbinom(2,1,0.5)) }
    else if (geno_mat_f0[p1,j]==1 && geno_mat_f0[p2,j]==2) { geno_mat_f1test[i,j] <- rbinom(1,1,0.5) + 1 }
    else if (geno_mat_f0[p1,j]==2 && geno_mat_f0[p2,j]==0) { geno_mat_f1test[i,j] <- 1 }
    else if (geno_mat_f0[p1,j]==2 && geno_mat_f0[p2,j]==1) { geno_mat_f1test[i,j] <- rbinom(1,1,0.5) + 1 }
    else if (geno_mat_f0[p1,j]==2 && geno_mat_f0[p2,j]==2) { geno_mat_f1test[i,j] <- 2 }
  }
}

# now create the genotype matrix for the test parents and offspring
# BUT ONLY THE PARENTS THAT WERE ACTUALLY CROSSED
crossed_parents <- unique(c(true_parents$parent1[1001:1095], true_parents$parent2[1001:1095]))
geno_mat_crossed_parents <- geno_mat_testparents[rownames(geno_mat_testparents) %in% crossed_parents, ]
dim(geno_mat_crossed_parents) # 68 parents

dim(geno_mat_f1test) # 95
dim(geno_mat_crossed_parents) # 68
geno_mat_bestcase <- rbind(geno_mat_crossed_parents, geno_mat_f1test)
dim(geno_mat_bestcase) # 163 inds. i mean we have no idea what the test group 
# SHOULD BE optimally... but let's push on with this in the morning.

# so now we have the best case set up. run sequoia on that. we have the true parents
# to compare it to. then we can do the 100% ptps with some extra (reminiscent of 
# our actual test group situation) if we use all of the 112 parents instead of 
# just the 68, or some different levels of that (make sure to calculate %extraparents),
# and we can start removing some of those 68 to get at <100% ptps

# first thing is to come up with a table that we can use to summarize the results
# of each situation, then we can let it rip.

##### ---- ---- #####

