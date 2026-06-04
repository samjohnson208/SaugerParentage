##### ---- sj_postdefense_simulations_trial1.R ---- #####
### by SPJ 060226 ###

library(sequoia)
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
nloci <- 25000
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
dim(geno_mat_f0f1) # 1000 inds, 25000 loci

##### ---- ---- #####

dim(geno_mat_f0f1)

table(table(unique(true_parents$parent1)))
table(table(unique(true_parents$parent2))) # 206 and 214 = 420
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
length(test_f0parents) # 112. good.

# create genotype matrix of broodstock contributors
geno_mat_testparents <- geno_mat_f0[rownames(geno_mat_f0) %in% test_f0parents, ]
dim(geno_mat_testparents) # 112 x 25000. good.

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

# now go into the true_parents (pedigree) dataframe, add the 95 test offspring rows
new_rows <- data.frame(
  offspring = rep(NA_character_, 95),
  parent1 = rep(NA_character_, 95),
  parent2 = rep(NA_character_, 95),
  sex = rep(NA_character_, 95),
  parent_cross = rep(NA_character_, 95)
)
# bind it with the 95 rows of NA
true_parents <- rbind(true_parents, new_rows)

# reset rownames
rownames(true_parents) <- NULL

# now paste in the pattern for test f1 offspring
true_parents$offspring[1001:1095] <- paste0("f1test_", 1:95)

# then for each one...
for(i in 1:95){
  row <- 1000 + i # define the individual that we're focusing on
  
  family <- family_draws[i] # grab that entry from the families vector (95 long)
  
  # then store the parents from the family in the true_parents dataframe
  true_parents$parent1[row] <- test_crosses$male[family]
  true_parents$parent2[row] <- test_crosses$female[family]
  
  # not sure why this has to be so complicated but it does... grab the parent
  # for the focal ind from true_parents.
  p1 <- true_parents$parent1[row]
  p2 <- true_parents$parent2[row]
  
  # and grab its number
  id1 <- as.numeric(sub("f0_", "", p1))
  id2 <- as.numeric(sub("f0_", "", p2))
  
  # then sort those numbers and paste them into parent_cross
  true_parents$parent_cross[row] <- paste(
    paste0("f0_", pmin(id1, id2)),
    paste0("f0_", pmax(id1, id2)),
    sep="__")
}

# now we assemble the genotype matrix for the test f1s
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
dim(geno_mat_crossed_parents) # 65 parents

dim(geno_mat_f1test) # 95
dim(geno_mat_crossed_parents) # 65
geno_mat_bestcase <- rbind(geno_mat_crossed_parents, geno_mat_f1test)
dim(geno_mat_bestcase) # 160 inds. i mean we have no idea what the test group 
# SHOULD BE optimally... but let's push on with this in the morning.

# so now we have the best case set up. run sequoia on that. we have the true parents
# to compare it to. then we can do the 100% ptps with some extra (reminiscent of 
# our actual test group situation) if we use all of the 112 parents instead of 
# just the 68, or some different levels of that (make sure to calculate %extraparents),
# and we can start removing some of those 68 to get at <100% ptps

# first thing is to come up with a table that we can use to summarize the results
# of each situation, then we can let it rip.

##### ---- ---- #####

######### ----- TESTING SEQUOIA PERFORMANCE ON SEVERAL DATASETS ----- #########

##### ---- create geno_mats ---- #####
##### ---- best case ---- #####

# filter the geno mat to maf 30
freqs_bestcase <- data.frame(site = 1:ncol(geno_mat_bestcase), site_mean = NA, site_mean_d2 = NA, fold_site_mean_d2 = NA)

for (i in 1:ncol(geno_mat_bestcase)) {
  freqs_bestcase$site_mean[i] <- mean(geno_mat_bestcase[, i], na.rm = TRUE)
}

for (i in 1:nrow(freqs_bestcase)) {
  m <- freqs_bestcase$site_mean[i]
  d2 <- m/2
  freqs_bestcase$site_mean_d2[i] <- d2
}

for (i in 1:nrow(freqs_bestcase)) {
  f <- freqs_bestcase$site_mean_d2[i]
  if (f > 0.5){
    freqs_bestcase$fold_site_mean_d2[i] <- 1-f
  } else {
    freqs_bestcase$fold_site_mean_d2[i] <- f
  }
}

sites_to_keep <- freqs_bestcase$site_mean_d2 >= 0.3 &
                 freqs_bestcase$site_mean_d2 <= 0.7

geno_mat_bestcase_filt <- geno_mat_bestcase[, sites_to_keep]
dim(geno_mat_bestcase_filt) # 160 inds and 1802 sites, that's so awesome.

# but we need 940, just sample the one time so all datasets have the same snps
# keep_cols <- sort(sample(ncol(geno_mat_bestcase_filt), 940))
geno_mat_bestcase_filt <- geno_mat_bestcase_filt[, keep_cols] 
dim(geno_mat_bestcase_filt)
# ^ PASTE JUST THIS WHEN DOING OTHER GROUPS
##### ---- 75% ptps ---- #####
dim(geno_mat_bestcase) # 65 parents, 95 offspring

# want to take it to 49 parents
dim(geno_mat_crossed_parents)
parents_to_keep <- sample(nrow(geno_mat_crossed_parents), 49, replace = FALSE)

geno_mat_75ptps <- geno_mat_bestcase[parents_to_keep, ]
geno_mat_75ptps <- rbind(geno_mat_75ptps, geno_mat_f1test)
dim(geno_mat_75ptps) # 144 = 95 + 49, okay, checks out.

# filter for maf
geno_mat_75ptps_filt <- geno_mat_75ptps[, keep_cols] 
dim(geno_mat_75ptps_filt) # 144 940. nice.

##### ---- 50% ptps ---- #####
dim(geno_mat_bestcase) # 65 parents, 95 offspring

# want to take it to 33 parents
dim(geno_mat_crossed_parents)
parents_to_keep <- sample(nrow(geno_mat_crossed_parents), 33, replace = FALSE)

geno_mat_50ptps <- geno_mat_bestcase[parents_to_keep, ]
geno_mat_50ptps <- rbind(geno_mat_50ptps, geno_mat_f1test)
dim(geno_mat_50ptps) # 128 = 95 + 33, okay, checks out.

# filter for maf
geno_mat_50ptps_filt <- geno_mat_50ptps[, keep_cols] 
dim(geno_mat_50ptps_filt) # 128 940. nice.
##### ---- 25% ptps ---- #####
dim(geno_mat_bestcase) # 65 parents, 95 offspring

# want to take it to 17 parents
dim(geno_mat_crossed_parents)
parents_to_keep <- sample(nrow(geno_mat_crossed_parents), 17, replace = FALSE)

geno_mat_25ptps <- geno_mat_bestcase[parents_to_keep, ]
geno_mat_25ptps <- rbind(geno_mat_25ptps, geno_mat_f1test)
dim(geno_mat_25ptps) # 112 = 95 + 17, okay, checks out.

# filter for maf
geno_mat_25ptps_filt <- geno_mat_25ptps[, keep_cols] 
dim(geno_mat_25ptps_filt) # 112 940. nice.
##### ---- 100% ptps + 25% extra ---- #####
# who was NOT crossed?
parents_not_crossed <- rownames(geno_mat_f0)[!(rownames(geno_mat_f0) %in% crossed_parents)]
# filter the f0 geno mat for those individuals
geno_mat_parents_not_crossed <- geno_mat_f0[rownames(geno_mat_f0) %in% parents_not_crossed,]

# sample 17 extra parents...
extra_parents <- geno_mat_parents_not_crossed[sample(nrow(geno_mat_parents_not_crossed), 17, replace = FALSE),]

# ... and add them to the best case scenario genomat
geno_mat_25extra <- rbind(extra_parents, geno_mat_bestcase)
dim(geno_mat_25extra) # 177 inds, 25000 snps. good.

# then filter that genomat for the snps we're interested in.
geno_mat_25extra_filt <- geno_mat_25extra[, keep_cols] 
dim(geno_mat_25extra_filt) # 177 and 940. good.

##### ---- 100% ptps + 50% extra ---- #####
# who was NOT crossed?
parents_not_crossed <- rownames(geno_mat_f0)[!(rownames(geno_mat_f0) %in% crossed_parents)]
# filter the f0 geno mat for those individuals
geno_mat_parents_not_crossed <- geno_mat_f0[rownames(geno_mat_f0) %in% parents_not_crossed,]

# sample 33 extra parents...
extra_parents <- geno_mat_parents_not_crossed[sample(nrow(geno_mat_parents_not_crossed), 33, replace = FALSE),]

# ... and add them to the best case scenario genomat
geno_mat_50extra <- rbind(extra_parents, geno_mat_bestcase)
dim(geno_mat_50extra) # 193 inds, 25000 snps. good.

# then filter that genomat for the snps we're interested in.
geno_mat_50extra_filt <- geno_mat_50extra[, keep_cols] 
dim(geno_mat_50extra_filt) # 193 and 940. good.
##### ---- 100% ptps + 75% extra ---- #####
# who was NOT crossed?
parents_not_crossed <- rownames(geno_mat_f0)[!(rownames(geno_mat_f0) %in% crossed_parents)]
# filter the f0 geno mat for those individuals
geno_mat_parents_not_crossed <- geno_mat_f0[rownames(geno_mat_f0) %in% parents_not_crossed,]

# sample 49 extra parents...
extra_parents <- geno_mat_parents_not_crossed[sample(nrow(geno_mat_parents_not_crossed), 49, replace = FALSE),]

# ... and add them to the best case scenario genomat
geno_mat_75extra <- rbind(extra_parents, geno_mat_bestcase)
dim(geno_mat_75extra) # 209 inds, 25000 snps. good.

# then filter that genomat for the snps we're interested in.
geno_mat_75extra_filt <- geno_mat_75extra[, keep_cols] 
dim(geno_mat_75extra_filt) # 209 and 940. good.
##### ---- 100% ptps + 100% extra ---- #####
# who was NOT crossed?
parents_not_crossed <- rownames(geno_mat_f0)[!(rownames(geno_mat_f0) %in% crossed_parents)]
# filter the f0 geno mat for those individuals
geno_mat_parents_not_crossed <- geno_mat_f0[rownames(geno_mat_f0) %in% parents_not_crossed,]

# sample 65 extra parents...
extra_parents <- geno_mat_parents_not_crossed[sample(nrow(geno_mat_parents_not_crossed), 65, replace = FALSE),]

# ... and add them to the best case scenario genomat
geno_mat_100extra <- rbind(extra_parents, geno_mat_bestcase)
dim(geno_mat_100extra) # 225 inds, 25000 snps. good.

# then filter that genomat for the snps we're interested in.
geno_mat_100extra_filt <- geno_mat_100extra[, keep_cols] 
dim(geno_mat_100extra_filt) # 225 and 940. good.
##### ---- ---- #####

#### ---- create LH dataframes ---- #####
##### ---- best case ---- #####
LH_bestcase <- data.frame(ID = rownames(geno_mat_bestcase),
                          Sex = NA,
                          BirthYear = NA)

LH_bestcase$BirthYear[1:65] <- 0
LH_bestcase$BirthYear[66:nrow(LH_bestcase)] <- 1

LH_bestcase$Sex <- true_parents$sex[match(LH_bestcase$ID, true_parents$offspring)]
LH_bestcase$Sex[LH_bestcase$Sex == "M"] <- 2
LH_bestcase$Sex[LH_bestcase$Sex == "F"] <- 1
LH_bestcase$Sex[66:nrow(LH_bestcase)] <- 3
##### ---- 75% ptps ---- #####
LH_75ptps <- data.frame(ID = rownames(geno_mat_75ptps),
                          Sex = NA,
                          BirthYear = NA)

LH_75ptps$BirthYear[1:49] <- 0
LH_75ptps$BirthYear[50:nrow(LH_75ptps)] <- 1

LH_75ptps$Sex <- true_parents$sex[match(LH_75ptps$ID, true_parents$offspring)]
LH_75ptps$Sex[LH_75ptps$Sex == "M"] <- 2
LH_75ptps$Sex[LH_75ptps$Sex == "F"] <- 1
LH_75ptps$Sex[50:nrow(LH_75ptps)] <- 3
##### ---- 50% ptps ---- #####
LH_50ptps <- data.frame(ID = rownames(geno_mat_50ptps),
                        Sex = NA,
                        BirthYear = NA)

LH_50ptps$BirthYear[1:33] <- 0
LH_50ptps$BirthYear[34:nrow(LH_50ptps)] <- 1

LH_50ptps$Sex <- true_parents$sex[match(LH_50ptps$ID, true_parents$offspring)]
LH_50ptps$Sex[LH_50ptps$Sex == "M"] <- 2
LH_50ptps$Sex[LH_50ptps$Sex == "F"] <- 1
LH_50ptps$Sex[34:nrow(LH_50ptps)] <- 3
##### ---- 25% ptps ---- #####
LH_25ptps <- data.frame(ID = rownames(geno_mat_25ptps),
                        Sex = NA,
                        BirthYear = NA)

LH_25ptps$BirthYear[1:17] <- 0
LH_25ptps$BirthYear[18:nrow(LH_25ptps)] <- 1

LH_25ptps$Sex <- true_parents$sex[match(LH_25ptps$ID, true_parents$offspring)]
LH_25ptps$Sex[LH_25ptps$Sex == "M"] <- 2
LH_25ptps$Sex[LH_25ptps$Sex == "F"] <- 1
LH_25ptps$Sex[18:nrow(LH_25ptps)] <- 3
##### ---- 100% ptps + 25% extra ---- #####
LH_25extra <- data.frame(ID = rownames(geno_mat_25extra_filt),
                          Sex = NA,
                          BirthYear = NA)

LH_25extra$BirthYear[1:82] <- 0
LH_25extra$BirthYear[83:nrow(LH_25extra)] <- 1

LH_25extra$Sex <- true_parents$sex[match(LH_25extra$ID, true_parents$offspring)]
LH_25extra$Sex[LH_25extra$Sex == "M"] <- 2
LH_25extra$Sex[LH_25extra$Sex == "F"] <- 1
LH_25extra$Sex[83:nrow(LH_25extra)] <- 3

##### ---- 100% ptps + 50% extra ---- #####
LH_50extra <- data.frame(ID = rownames(geno_mat_50extra_filt),
                         Sex = NA,
                         BirthYear = NA)

LH_50extra$BirthYear[1:98] <- 0
LH_50extra$BirthYear[99:nrow(LH_50extra)] <- 1

LH_50extra$Sex <- true_parents$sex[match(LH_50extra$ID, true_parents$offspring)]
LH_50extra$Sex[LH_50extra$Sex == "M"] <- 2
LH_50extra$Sex[LH_50extra$Sex == "F"] <- 1
LH_50extra$Sex[99:nrow(LH_50extra)] <- 3
##### ---- 100% ptps + 75% extra ---- #####
LH_75extra <- data.frame(ID = rownames(geno_mat_75extra_filt),
                         Sex = NA,
                         BirthYear = NA)

LH_75extra$BirthYear[1:114] <- 0
LH_75extra$BirthYear[115:nrow(LH_75extra)] <- 1

LH_75extra$Sex <- true_parents$sex[match(LH_75extra$ID, true_parents$offspring)]
LH_75extra$Sex[LH_75extra$Sex == "M"] <- 2
LH_75extra$Sex[LH_75extra$Sex == "F"] <- 1
LH_75extra$Sex[115:nrow(LH_75extra)] <- 3
##### ---- 100% ptps + 100% extra ---- #####
LH_100extra <- data.frame(ID = rownames(geno_mat_100extra_filt),
                         Sex = NA,
                         BirthYear = NA)

LH_100extra$BirthYear[1:130] <- 0
LH_100extra$BirthYear[131:nrow(LH_100extra)] <- 1

LH_100extra$Sex <- true_parents$sex[match(LH_100extra$ID, true_parents$offspring)]
LH_100extra$Sex[LH_100extra$Sex == "M"] <- 2
LH_100extra$Sex[LH_100extra$Sex == "F"] <- 1
LH_100extra$Sex[131:nrow(LH_100extra)] <- 3
##### ---- ---- #####

##### ---- run sequoia(), GetMaybeRel() ---- #####
# all params the same as for empirical data, but error is commented out.
##### ---- best case ---- #####
seq_bestcase <- sequoia(GenoM = geno_mat_bestcase_filt,
                  LifeHistData = LH_bestcase,
                  Module = "ped",
                  # Err = errM,
                  Complex = "full",
                  Herm = "no",
                  UseAge = "yes",
                  args.AP=list(Discrete = TRUE, 
                               MinAgeParent = 1, MaxAgeParent = 1),
                  CalcLLR = TRUE,
                  StrictGenoCheck = TRUE,
                  DummyPrefix = c("F", "M"),
                  Tfilter = -2,
                  Tassign = 0.5)
# assigned 95 dams and 95 sires to 160 + 0 individuals (real + dummy)

gmr_bestcase <- GetMaybeRel(GenoM = geno_mat_bestcase_filt,
                      SeqList = seq_bestcase,
                      AgePrior = seq_bestcase[["AgePriors"]],
                      # Err = errM,
                      Module = "ped",
                      Complex = "full",
                      LifeHistData = LH_bestcase,
                      Herm = "no",
                      quiet = FALSE,
                      Tfilter = -2,
                      Tassign = 0.5,
                      MaxPairs = 7 * nrow(geno_mat_bestcase_filt))
# Found 0 likely parent-offspring pairs, and 0, other non-assigned pairs of possible relatives
##### ---- 75% ptps---- #####
seq_75ptps <- sequoia(GenoM = geno_mat_75ptps_filt,
                        LifeHistData = LH_75ptps,
                        Module = "ped",
                        # Err = errM,
                        Complex = "full",
                        Herm = "no",
                        UseAge = "yes",
                        args.AP=list(Discrete = TRUE, 
                                     MinAgeParent = 1, MaxAgeParent = 1),
                        CalcLLR = TRUE,
                        StrictGenoCheck = TRUE,
                        DummyPrefix = c("F", "M"),
                        Tfilter = -2,
                        Tassign = 0.5)
# assigned 89 dams and 90 sires to 144 + 11 individuals

gmr_75ptps <- GetMaybeRel(GenoM = geno_mat_75ptps_filt,
                            SeqList = seq_75ptps,
                            AgePrior = seq_75ptps[["AgePriors"]],
                            # Err = errM,
                            Module = "ped",
                            Complex = "full",
                            LifeHistData = LH_75ptps,
                            Herm = "no",
                            quiet = FALSE,
                            Tfilter = -2,
                            Tassign = 0.5,
                            MaxPairs = 7 * nrow(geno_mat_75ptps_filt))
# Found 0 likely parent-offspring pairs, and 5, other non-assigned pairs of possible relatives
##### ---- 50% ptps---- #####
seq_50ptps <- sequoia(GenoM = geno_mat_50ptps_filt,
                      LifeHistData = LH_50ptps,
                      Module = "ped",
                      # Err = errM,
                      Complex = "full",
                      Herm = "no",
                      UseAge = "yes",
                      args.AP=list(Discrete = TRUE, 
                                   MinAgeParent = 1, MaxAgeParent = 1),
                      CalcLLR = TRUE,
                      StrictGenoCheck = TRUE,
                      DummyPrefix = c("F", "M"),
                      Tfilter = -2,
                      Tassign = 0.5)
# assigned 89 dams and 85 sires to 128 + 21 individuals (real + dummy)

gmr_50ptps <- GetMaybeRel(GenoM = geno_mat_50ptps_filt,
                          SeqList = seq_50ptps,
                          AgePrior = seq_50ptps[["AgePriors"]],
                          # Err = errM,
                          Module = "ped",
                          Complex = "full",
                          LifeHistData = LH_50ptps,
                          Herm = "no",
                          quiet = FALSE,
                          Tfilter = -2,
                          Tassign = 0.5,
                          MaxPairs = 7 * nrow(geno_mat_50ptps_filt))
# Found 0 likely parent-offspring pairs, and 3, other non-assigned pairs of possible relatives
##### ---- 25% ptps---- ##### 
seq_25ptps <- sequoia(GenoM = geno_mat_25ptps_filt,
                      LifeHistData = LH_25ptps,
                      Module = "ped",
                      # Err = errM,
                      Complex = "full",
                      Herm = "no",
                      UseAge = "yes",
                      args.AP=list(Discrete = TRUE, 
                                   MinAgeParent = 1, MaxAgeParent = 1),
                      CalcLLR = TRUE,
                      StrictGenoCheck = TRUE,
                      DummyPrefix = c("F", "M"),
                      Tfilter = -2,
                      Tassign = 0.5)
# assigned 81 dams and 80 sires to 112 + 34 individuals (real + dummy)

gmr_25ptps <- GetMaybeRel(GenoM = geno_mat_25ptps_filt,
                          SeqList = seq_25ptps,
                          AgePrior = seq_25ptps[["AgePriors"]],
                          # Err = errM,
                          Module = "ped",
                          Complex = "full",
                          LifeHistData = LH_25ptps,
                          Herm = "no",
                          quiet = FALSE,
                          Tfilter = -2,
                          Tassign = 0.5,
                          MaxPairs = 7 * nrow(geno_mat_25ptps_filt))
# Found 0 likely parent-offspring pairs, and 26, other non-assigned pairs of possible relatives
##### ---- 100% + 25% extra---- ##### 
seq_25extra <- sequoia(GenoM = geno_mat_25extra_filt,
                      LifeHistData = LH_25extra,
                      Module = "ped",
                      # Err = errM,
                      Complex = "full",
                      Herm = "no",
                      UseAge = "yes",
                      args.AP=list(Discrete = TRUE, 
                                   MinAgeParent = 1, MaxAgeParent = 1),
                      CalcLLR = TRUE,
                      StrictGenoCheck = TRUE,
                      DummyPrefix = c("F", "M"),
                      Tfilter = -2,
                      Tassign = 0.5)
# assigned 95 dams and 95 sires to 177 + 0 individuals (real + dummy) 

gmr_25extra <- GetMaybeRel(GenoM = geno_mat_25extra_filt,
                          SeqList = seq_25extra,
                          AgePrior = seq_25extra[["AgePriors"]],
                          # Err = errM,
                          Module = "ped",
                          Complex = "full",
                          LifeHistData = LH_25extra,
                          Herm = "no",
                          quiet = FALSE,
                          Tfilter = -2,
                          Tassign = 0.5,
                          MaxPairs = 7 * nrow(geno_mat_25extra_filt))
# Found 0 likely parent-offspring pairs, and 1, other non-assigned pairs of possible relatives
##### ---- 100% + 50% extra---- ##### 
seq_50extra <- sequoia(GenoM = geno_mat_50extra_filt,
                       LifeHistData = LH_50extra,
                       Module = "ped",
                       # Err = errM,
                       Complex = "full",
                       Herm = "no",
                       UseAge = "yes",
                       args.AP=list(Discrete = TRUE, 
                                    MinAgeParent = 1, MaxAgeParent = 1),
                       CalcLLR = TRUE,
                       StrictGenoCheck = TRUE,
                       DummyPrefix = c("F", "M"),
                       Tfilter = -2,
                       Tassign = 0.5)
# assigned 92 dams and 92 sires to 193 + 4 individuals (real + dummy) 

gmr_50extra <- GetMaybeRel(GenoM = geno_mat_50extra_filt,
                           SeqList = seq_50extra,
                           AgePrior = seq_50extra[["AgePriors"]],
                           # Err = errM,
                           Module = "ped",
                           Complex = "full",
                           LifeHistData = LH_50extra,
                           Herm = "no",
                           quiet = FALSE,
                           Tfilter = -2,
                           Tassign = 0.5,
                           MaxPairs = 7 * nrow(geno_mat_50extra_filt))
# Found 6 likely parent-offspring pairs, and 19, other non-assigned pairs of possible relatives
# Found 3 parent-parent-offspring trios
##### ---- 100% + 75% extra---- ##### 
seq_75extra <- sequoia(GenoM = geno_mat_75extra_filt,
                       LifeHistData = LH_75extra,
                       Module = "ped",
                       # Err = errM,
                       Complex = "full",
                       Herm = "no",
                       UseAge = "yes",
                       args.AP=list(Discrete = TRUE, 
                                    MinAgeParent = 1, MaxAgeParent = 1),
                       CalcLLR = TRUE,
                       StrictGenoCheck = TRUE,
                       DummyPrefix = c("F", "M"),
                       Tfilter = -2,
                       Tassign = 0.5)
# assigned 95 dams and 95 sires to 209 + 0 individuals (real + dummy) 

gmr_75extra <- GetMaybeRel(GenoM = geno_mat_75extra_filt,
                           SeqList = seq_75extra,
                           AgePrior = seq_75extra[["AgePriors"]],
                           # Err = errM,
                           Module = "ped",
                           Complex = "full",
                           LifeHistData = LH_75extra,
                           Herm = "no",
                           quiet = FALSE,
                           Tfilter = -2,
                           Tassign = 0.5,
                           MaxPairs = 7 * nrow(geno_mat_75extra_filt))
# Found 0 likely parent-offspring pairs, and 3, other non-assigned pairs of possible relatives
##### ---- 100% + 100% extra---- ##### 
seq_100extra <- sequoia(GenoM = geno_mat_100extra_filt,
                       LifeHistData = LH_100extra,
                       Module = "ped",
                       # Err = errM,
                       Complex = "full",
                       Herm = "no",
                       UseAge = "yes",
                       args.AP=list(Discrete = TRUE, 
                                    MinAgeParent = 1, MaxAgeParent = 1),
                       CalcLLR = TRUE,
                       StrictGenoCheck = TRUE,
                       DummyPrefix = c("F", "M"),
                       Tfilter = -2,
                       Tassign = 0.5)
#  assigned 95 dams and 95 sires to 225 + 0 individuals (real + dummy) 

gmr_100extra <- GetMaybeRel(GenoM = geno_mat_100extra_filt,
                           SeqList = seq_100extra,
                           AgePrior = seq_100extra[["AgePriors"]],
                           # Err = errM,
                           Module = "ped",
                           Complex = "full",
                           LifeHistData = LH_100extra,
                           Herm = "no",
                           quiet = FALSE,
                           Tfilter = -2,
                           Tassign = 0.5,
                           MaxPairs = 7 * nrow(geno_mat_100extra_filt))
# Found 0 likely parent-offspring pairs, and 5, other non-assigned pairs of possible relatives
##### ---- ---- #####

##### ---- Pairs_bestcase ---- #####
IDs <- rownames(geno_mat_bestcase)
Pairs_bestcase <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_bestcase <- Pairs_bestcase %>% 
  left_join(LH_bestcase %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_bestcase %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_bestcase$AgeDif <- Pairs_bestcase$BY2 - Pairs_bestcase$BY1

Pairs_bestcase$focal <- "U"

dim(Pairs_bestcase)
##### ---- Pairs_75ptps ---- #####
IDs <- rownames(geno_mat_75ptps_filt)
Pairs_75ptps <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_75ptps <- Pairs_75ptps %>% 
  left_join(LH_75ptps %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_75ptps %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_75ptps$AgeDif <- Pairs_75ptps$BY2 - Pairs_75ptps$BY1

Pairs_75ptps$focal <- "U"

dim(Pairs_75ptps)

##### ---- Pairs_50ptps ---- #####
IDs <- rownames(geno_mat_50ptps_filt)
Pairs_50ptps <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_50ptps <- Pairs_50ptps %>% 
  left_join(LH_50ptps %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_50ptps %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_50ptps$AgeDif <- Pairs_50ptps$BY2 - Pairs_50ptps$BY1

Pairs_50ptps$focal <- "U"

dim(Pairs_50ptps)
##### ---- Pairs_25ptps ---- #####
IDs <- rownames(geno_mat_25ptps_filt)
Pairs_25ptps <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_25ptps <- Pairs_25ptps %>% 
  left_join(LH_25ptps %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_25ptps %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_25ptps$AgeDif <- Pairs_25ptps$BY2 - Pairs_25ptps$BY1

Pairs_25ptps$focal <- "U"

dim(Pairs_25ptps)
##### ---- Pairs_25extra ---- #####
IDs <- rownames(geno_mat_25extra_filt)
Pairs_25extra <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_25extra <- Pairs_25extra %>% 
  left_join(LH_25extra %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_25extra %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_25extra$AgeDif <- Pairs_25extra$BY2 - Pairs_25extra$BY1

Pairs_25extra$focal <- "U"

dim(Pairs_25extra)
##### ---- Pairs_50extra ---- #####
IDs <- rownames(geno_mat_50extra_filt)
Pairs_50extra <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_50extra <- Pairs_50extra %>% 
  left_join(LH_50extra %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_50extra %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_50extra$AgeDif <- Pairs_50extra$BY2 - Pairs_50extra$BY1

Pairs_50extra$focal <- "U"

dim(Pairs_50extra)
##### ---- Pairs_75extra ---- #####
IDs <- rownames(geno_mat_75extra_filt)
Pairs_75extra <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_75extra <- Pairs_75extra %>% 
  left_join(LH_75extra %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_75extra %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_75extra$AgeDif <- Pairs_75extra$BY2 - Pairs_75extra$BY1

Pairs_75extra$focal <- "U"

dim(Pairs_75extra)
##### ---- Pairs_100extra ---- #####
IDs <- rownames(geno_mat_100extra_filt)
Pairs_100extra <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_100extra <- Pairs_100extra %>% 
  left_join(LH_100extra %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_100extra %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_100extra$AgeDif <- Pairs_100extra$BY2 - Pairs_100extra$BY1

Pairs_100extra$focal <- "U"

dim(Pairs_100extra)
##### ---- ---- #####

##### ---- Getting LLRs and probs for all relationships ---- #####
##### ---- PairLL_bestcase -> prob_pairs_bestcase ---- #####
PairLL_bestcase <- CalcPairLL(Pairs = Pairs_bestcase,
                          GenoM = geno_mat_bestcase_filt,
                          LifeHistData = LH_bestcase,
                          AgePrior = seq_bestcase[["AgePriors"]],
                          Module = "ped",
                          Complex = "full",
                          Herm = 'no',
                          InclDup = FALSE,
                          # Err = errM,
                          Tassign = 0.5,
                          Tfilter = -2,
                          quiet = FALSE,
                          Plot = TRUE)
prob_pairs_bestcase <- plyr::aaply(as.matrix(PairLL_bestcase[,10:16]), .margin = 1, LLtoProb)
prob_pairs_bestcase <- cbind(PairLL_bestcase[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_bestcase)
##### ---- PairLL_75ptps -> prob_pairs_75ptps ---- #####
PairLL_75ptps <- CalcPairLL(Pairs = Pairs_75ptps,
                              GenoM = geno_mat_75ptps_filt,
                              LifeHistData = LH_75ptps,
                              AgePrior = seq_75ptps[["AgePriors"]],
                              Module = "ped",
                              Complex = "full",
                              Herm = 'no',
                              InclDup = FALSE,
                              # Err = errM,
                              Tassign = 0.5,
                              Tfilter = -2,
                              quiet = FALSE,
                              Plot = TRUE)
prob_pairs_75ptps <- plyr::aaply(as.matrix(PairLL_75ptps[,10:16]), .margin = 1, LLtoProb)
prob_pairs_75ptps <- cbind(PairLL_75ptps[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_75ptps)
##### ---- PairLL_50ptps -> prob_pairs_50ptps ---- #####
PairLL_50ptps <- CalcPairLL(Pairs = Pairs_50ptps,
                            GenoM = geno_mat_50ptps_filt,
                            LifeHistData = LH_50ptps,
                            AgePrior = seq_50ptps[["AgePriors"]],
                            Module = "ped",
                            Complex = "full",
                            Herm = 'no',
                            InclDup = FALSE,
                            # Err = errM,
                            Tassign = 0.5,
                            Tfilter = -2,
                            quiet = FALSE,
                            Plot = TRUE)
prob_pairs_50ptps <- plyr::aaply(as.matrix(PairLL_50ptps[,10:16]), .margin = 1, LLtoProb)
prob_pairs_50ptps <- cbind(PairLL_50ptps[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_50ptps)
##### ---- PairLL_25ptps -> prob_pairs_25ptps ---- #####
PairLL_25ptps <- CalcPairLL(Pairs = Pairs_25ptps,
                            GenoM = geno_mat_25ptps_filt,
                            LifeHistData = LH_25ptps,
                            AgePrior = seq_25ptps[["AgePriors"]],
                            Module = "ped",
                            Complex = "full",
                            Herm = 'no',
                            InclDup = FALSE,
                            # Err = errM,
                            Tassign = 0.5,
                            Tfilter = -2,
                            quiet = FALSE,
                            Plot = TRUE)
prob_pairs_25ptps <- plyr::aaply(as.matrix(PairLL_25ptps[,10:16]), .margin = 1, LLtoProb)
prob_pairs_25ptps <- cbind(PairLL_25ptps[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_25ptps)
##### ---- PairLL_25extra -> prob_pairs_25extra ---- #####
PairLL_25extra <- CalcPairLL(Pairs = Pairs_25extra,
                            GenoM = geno_mat_25extra_filt,
                            LifeHistData = LH_25extra,
                            AgePrior = seq_25extra[["AgePriors"]],
                            Module = "ped",
                            Complex = "full",
                            Herm = 'no',
                            InclDup = FALSE,
                            # Err = errM,
                            Tassign = 0.5,
                            Tfilter = -2,
                            quiet = FALSE,
                            Plot = TRUE)
prob_pairs_25extra <- plyr::aaply(as.matrix(PairLL_25extra[,10:16]), .margin = 1, LLtoProb)
prob_pairs_25extra <- cbind(PairLL_25extra[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_25extra)
##### ---- PairLL_50extra -> prob_pairs_50extra ---- #####
PairLL_50extra <- CalcPairLL(Pairs = Pairs_50extra,
                             GenoM = geno_mat_50extra_filt,
                             LifeHistData = LH_50extra,
                             AgePrior = seq_50extra[["AgePriors"]],
                             Module = "ped",
                             Complex = "full",
                             Herm = 'no',
                             InclDup = FALSE,
                             # Err = errM,
                             Tassign = 0.5,
                             Tfilter = -2,
                             quiet = FALSE,
                             Plot = TRUE)
prob_pairs_50extra <- plyr::aaply(as.matrix(PairLL_50extra[,10:16]), .margin = 1, LLtoProb)
prob_pairs_50extra <- cbind(PairLL_50extra[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_50extra)
##### ---- PairLL_75extra -> prob_pairs_75extra ---- #####
PairLL_75extra <- CalcPairLL(Pairs = Pairs_75extra,
                             GenoM = geno_mat_75extra_filt,
                             LifeHistData = LH_75extra,
                             AgePrior = seq_75extra[["AgePriors"]],
                             Module = "ped",
                             Complex = "full",
                             Herm = 'no',
                             InclDup = FALSE,
                             # Err = errM,
                             Tassign = 0.5,
                             Tfilter = -2,
                             quiet = FALSE,
                             Plot = TRUE)
prob_pairs_75extra <- plyr::aaply(as.matrix(PairLL_75extra[,10:16]), .margin = 1, LLtoProb)
prob_pairs_75extra <- cbind(PairLL_75extra[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_75extra)
##### ---- PairLL_100extra -> prob_pairs_100extra ---- #####
PairLL_100extra <- CalcPairLL(Pairs = Pairs_100extra,
                             GenoM = geno_mat_100extra_filt,
                             LifeHistData = LH_100extra,
                             AgePrior = seq_100extra[["AgePriors"]],
                             Module = "ped",
                             Complex = "full",
                             Herm = 'no',
                             InclDup = FALSE,
                             # Err = errM,
                             Tassign = 0.5,
                             Tfilter = -2,
                             quiet = FALSE,
                             Plot = TRUE)
prob_pairs_100extra <- plyr::aaply(as.matrix(PairLL_100extra[,10:16]), .margin = 1, LLtoProb)
prob_pairs_100extra <- cbind(PairLL_100extra[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_100extra)
##### ---- ---- #####

##### ---- after creating those... ---- #####

##### ---- keep only distinct pairs ---- #####
prob_pairs_bestcase_unique <- prob_pairs_bestcase %>%
  mutate(
    id1_num = as.numeric(sub(".*_", "", ID1)),
    id2_num = as.numeric(sub(".*_", "", ID2)),
    Pair = ifelse(
      id1_num < id2_num,
      paste(ID1, ID2, sep = "__"),
      paste(ID2, ID1, sep = "__")
    )
  ) %>%
  distinct(Pair, .keep_all = TRUE) %>%
  select(-id1_num, -id2_num)

prob_pairs_75ptps_unique <- prob_pairs_75ptps %>%
  mutate(
    id1_num = as.numeric(sub(".*_", "", ID1)),
    id2_num = as.numeric(sub(".*_", "", ID2)),
    Pair = ifelse(
      id1_num < id2_num,
      paste(ID1, ID2, sep = "__"),
      paste(ID2, ID1, sep = "__")
    )
  ) %>%
  distinct(Pair, .keep_all = TRUE) %>%
  select(-id1_num, -id2_num)

prob_pairs_50ptps_unique <- prob_pairs_50ptps %>%
  mutate(
    id1_num = as.numeric(sub(".*_", "", ID1)),
    id2_num = as.numeric(sub(".*_", "", ID2)),
    Pair = ifelse(
      id1_num < id2_num,
      paste(ID1, ID2, sep = "__"),
      paste(ID2, ID1, sep = "__")
    )
  ) %>%
  distinct(Pair, .keep_all = TRUE) %>%
  select(-id1_num, -id2_num)

prob_pairs_25ptps_unique <- prob_pairs_25ptps %>%
  mutate(
    id1_num = as.numeric(sub(".*_", "", ID1)),
    id2_num = as.numeric(sub(".*_", "", ID2)),
    Pair = ifelse(
      id1_num < id2_num,
      paste(ID1, ID2, sep = "__"),
      paste(ID2, ID1, sep = "__")
    )
  ) %>%
  distinct(Pair, .keep_all = TRUE) %>%
  select(-id1_num, -id2_num)

prob_pairs_25extra_unique <- prob_pairs_25extra %>%
  mutate(
    id1_num = as.numeric(sub(".*_", "", ID1)),
    id2_num = as.numeric(sub(".*_", "", ID2)),
    Pair = ifelse(
      id1_num < id2_num,
      paste(ID1, ID2, sep = "__"),
      paste(ID2, ID1, sep = "__")
    )
  ) %>%
  distinct(Pair, .keep_all = TRUE) %>%
  select(-id1_num, -id2_num)

prob_pairs_50extra_unique <- prob_pairs_50extra %>%
  mutate(
    id1_num = as.numeric(sub(".*_", "", ID1)),
    id2_num = as.numeric(sub(".*_", "", ID2)),
    Pair = ifelse(
      id1_num < id2_num,
      paste(ID1, ID2, sep = "__"),
      paste(ID2, ID1, sep = "__")
    )
  ) %>%
  distinct(Pair, .keep_all = TRUE) %>%
  select(-id1_num, -id2_num)

prob_pairs_75extra_unique <- prob_pairs_75extra %>%
  mutate(
    id1_num = as.numeric(sub(".*_", "", ID1)),
    id2_num = as.numeric(sub(".*_", "", ID2)),
    Pair = ifelse(
      id1_num < id2_num,
      paste(ID1, ID2, sep = "__"),
      paste(ID2, ID1, sep = "__")
    )
  ) %>%
  distinct(Pair, .keep_all = TRUE) %>%
  select(-id1_num, -id2_num)

prob_pairs_100extra_unique <- prob_pairs_100extra %>%
  mutate(
    id1_num = as.numeric(sub(".*_", "", ID1)),
    id2_num = as.numeric(sub(".*_", "", ID2)),
    Pair = ifelse(
      id1_num < id2_num,
      paste(ID1, ID2, sep = "__"),
      paste(ID2, ID1, sep = "__")
    )
  ) %>%
  distinct(Pair, .keep_all = TRUE) %>%
  select(-id1_num, -id2_num)

##### ---- ---- #####

##### ---- assign generations to the pairs ---- #####
assign_gen <- function(x) {
  case_when(
    x %in% rownames(geno_mat_f0) ~ "F0",
    x %in% rownames(geno_mat_f1test) ~ "F1_Test",
    x %in% rownames(geno_mat_f1) ~ "F1",
    TRUE ~ NA_character_ # shouldn't be the case, but put an NA if an ind isn't there
  )
}

prob_pairs_bestcase_unique_gen <- prob_pairs_bestcase_unique %>%
  mutate(
    group_ind1 = assign_gen(ID1),
    group_ind2 = assign_gen(ID2),
    group_pair = paste(
      pmin(group_ind1, group_ind2),
      pmax(group_ind1, group_ind2),
      sep = "_"
    )
  )

prob_pairs_75ptps_unique_gen <- prob_pairs_75ptps_unique %>%
  mutate(
    group_ind1 = assign_gen(ID1),
    group_ind2 = assign_gen(ID2),
    group_pair = paste(
      pmin(group_ind1, group_ind2),
      pmax(group_ind1, group_ind2),
      sep = "_"
    )
  )

prob_pairs_50ptps_unique_gen <- prob_pairs_50ptps_unique %>%
  mutate(
    group_ind1 = assign_gen(ID1),
    group_ind2 = assign_gen(ID2),
    group_pair = paste(
      pmin(group_ind1, group_ind2),
      pmax(group_ind1, group_ind2),
      sep = "_"
    )
  )

prob_pairs_25ptps_unique_gen <- prob_pairs_25ptps_unique %>%
  mutate(
    group_ind1 = assign_gen(ID1),
    group_ind2 = assign_gen(ID2),
    group_pair = paste(
      pmin(group_ind1, group_ind2),
      pmax(group_ind1, group_ind2),
      sep = "_"
    )
  )

prob_pairs_25extra_unique_gen <- prob_pairs_25extra_unique %>%
  mutate(
    group_ind1 = assign_gen(ID1),
    group_ind2 = assign_gen(ID2),
    group_pair = paste(
      pmin(group_ind1, group_ind2),
      pmax(group_ind1, group_ind2),
      sep = "_"
    )
  )

prob_pairs_50extra_unique_gen <- prob_pairs_50extra_unique %>%
  mutate(
    group_ind1 = assign_gen(ID1),
    group_ind2 = assign_gen(ID2),
    group_pair = paste(
      pmin(group_ind1, group_ind2),
      pmax(group_ind1, group_ind2),
      sep = "_"
    )
  )

prob_pairs_75extra_unique_gen <- prob_pairs_75extra_unique %>%
  mutate(
    group_ind1 = assign_gen(ID1),
    group_ind2 = assign_gen(ID2),
    group_pair = paste(
      pmin(group_ind1, group_ind2),
      pmax(group_ind1, group_ind2),
      sep = "_"
    )
  )

prob_pairs_100extra_unique_gen <- prob_pairs_100extra_unique %>%
  mutate(
    group_ind1 = assign_gen(ID1),
    group_ind2 = assign_gen(ID2),
    group_pair = paste(
      pmin(group_ind1, group_ind2),
      pmax(group_ind1, group_ind2),
      sep = "_"
    )
  )


##### ---- ---- #####

##### ---- take just the POs ---- #####
PO_bestcase <- prob_pairs_bestcase_unique_gen %>% 
  filter(TopRel == "PO") %>% 
  filter(AgeDif == 1)

PO_75ptps <- prob_pairs_75ptps_unique_gen %>% 
  filter(TopRel == "PO") %>% 
  filter(AgeDif == 1)

PO_50ptps <- prob_pairs_50ptps_unique_gen %>% 
  filter(TopRel == "PO") %>% 
  filter(AgeDif == 1)

PO_25ptps <- prob_pairs_25ptps_unique_gen %>% 
  filter(TopRel == "PO") %>% 
  filter(AgeDif == 1)

PO_25extra <- prob_pairs_25extra_unique_gen %>% 
  filter(TopRel == "PO") %>% 
  filter(AgeDif == 1)

PO_50extra <- prob_pairs_50extra_unique_gen %>% 
  filter(TopRel == "PO") %>% 
  filter(AgeDif == 1)

PO_75extra <- prob_pairs_75extra_unique_gen %>% 
  filter(TopRel == "PO") %>% 
  filter(AgeDif == 1)

PO_100extra <- prob_pairs_100extra_unique_gen %>% 
  filter(TopRel == "PO") %>% 
  filter(AgeDif == 1)

##### ---- ---- #####

#### ---- how many assignments and accurate assignments in each situation? ---- ####

# here are the objects we're summarizing
PO_bestcase
PO_75ptps
PO_50ptps
PO_25ptps
PO_50extra
PO_75extra
PO_100extra

##### ---- PO_bestcase ---- #####
dim(PO_bestcase) # 190 assignments made
length(unique(PO_bestcase$ID2)) # 95 test f1's assigned
table(table(PO_bestcase$ID2)) # 95 assigned to two
length(unique(PO_bestcase$ID1)) # all 65 parents assigned
summary(PO_bestcase$PO) # all ones

PO_bestcase_counts <- PO_bestcase %>% 
  group_by(ID2) %>%
  mutate(n_parents = n()) %>%
  ungroup()

table(PO_bestcase_counts$n_parents) 
PO_bestcase_1 <- PO_bestcase_counts %>% filter(n_parents == 1) # 0
PO_bestcase_2 <- PO_bestcase_counts %>% filter(n_parents == 2) # 190

##### ---- PO_75ptps ---- #####
dim(PO_75ptps) # 136 assignments made
length(unique(PO_75ptps$ID2)) # 92 test f1's assigned
table(table(PO_75ptps$ID2)) # 44 assigned to two, 48 assigned to one
length(unique(PO_75ptps$ID1)) # all 49 parents assigned
summary(PO_75ptps$PO) # very high confidence

PO_75ptps_counts <- PO_75ptps %>% 
  group_by(ID2) %>%
  mutate(n_parents = n()) %>%
  ungroup()

table(PO_75ptps_counts$n_parents) 
PO_75ptps_1 <- PO_75ptps_counts %>% filter(n_parents == 1) # 48
PO_75ptps_2 <- PO_75ptps_counts %>% filter(n_parents == 2) # 44
##### ---- PO_50ptps ---- #####
dim(PO_50ptps) # 136 assignments made
length(unique(PO_50ptps$ID2)) # 72 test f1's assigned
table(table(PO_50ptps$ID2)) # 19 assigned to two, 53 assigned to one
length(unique(PO_50ptps$ID1)) # all 33 parents assigned
summary(PO_50ptps$PO) # very high confidence

PO_50ptps_counts <- PO_50ptps %>% 
  group_by(ID2) %>%
  mutate(n_parents = n()) %>%
  ungroup()

table(PO_50ptps_counts$n_parents) 
PO_50ptps_1 <- PO_50ptps_counts %>% filter(n_parents == 1) # 53
PO_50ptps_2 <- PO_50ptps_counts %>% filter(n_parents == 2) # 19
##### ---- PO_25ptps ---- #####
dim(PO_25ptps) # 55 assignments made
length(unique(PO_25ptps$ID2)) # 50 test f1's assigned
table(table(PO_25ptps$ID2)) # 5 assigned to two, 45 assigned to one
length(unique(PO_25ptps$ID1)) # all 17 parents assigned
summary(PO_25ptps$PO) # very high confidence

PO_25ptps_counts <- PO_25ptps %>% 
  group_by(ID2) %>%
  mutate(n_parents = n()) %>%
  ungroup()

table(PO_25ptps_counts$n_parents) 
PO_25ptps_1 <- PO_25ptps_counts %>% filter(n_parents == 1) # 45
PO_25ptps_2 <- PO_25ptps_counts %>% filter(n_parents == 2) # 5
##### ---- PO_25extra ---- #####
dim(PO_25extra) # 190 assignments made
length(unique(PO_25extra$ID2)) # 95 test f1's assigned
table(table(PO_25extra$ID2)) # 95 assigned to two
length(unique(PO_25extra$ID1)) # 65 of 65 parents assigned
summary(PO_25extra$PO) # very high confidence

PO_25extra_counts <- PO_25extra %>% 
  group_by(ID2) %>%
  mutate(n_parents = n()) %>%
  ungroup()

table(PO_25extra_counts$n_parents) 
PO_25extra_1 <- PO_25extra_counts %>% filter(n_parents == 1) # 45
PO_25extra_2 <- PO_25extra_counts %>% filter(n_parents == 2) # 5
##### ---- PO_50extra ---- #####
dim(PO_50extra) # 190 assignments made
length(unique(PO_50extra$ID2)) # 95 test f1's assigned
table(table(PO_50extra$ID2)) # 95 assigned to two
length(unique(PO_50extra$ID1)) # 65 of 65 parents assigned
summary(PO_50extra$PO) # very high confidence

PO_50extra_counts <- PO_50extra %>% 
  group_by(ID2) %>%
  mutate(n_parents = n()) %>%
  ungroup()

table(PO_50extra_counts$n_parents) 
PO_50extra_1 <- PO_50extra_counts %>% filter(n_parents == 1) # 0
PO_50extra_2 <- PO_50extra_counts %>% filter(n_parents == 2) # 95
##### ---- PO_75extra ---- #####
dim(PO_75extra) # 190 assignments made
length(unique(PO_75extra$ID2)) # 95 test f1's assigned
table(table(PO_75extra$ID2)) # 95 assigned to two
length(unique(PO_75extra$ID1)) # 65 of 65 parents assigned
summary(PO_75extra$PO) # very high confidence

PO_75extra_counts <- PO_75extra %>% 
  group_by(ID2) %>%
  mutate(n_parents = n()) %>%
  ungroup()

table(PO_75extra_counts$n_parents) 
PO_75extra_1 <- PO_75extra_counts %>% filter(n_parents == 1) # 0
PO_75extra_2 <- PO_75extra_counts %>% filter(n_parents == 2) # 95
##### ---- PO_100extra ---- #####
dim(PO_100extra) # 190 assignments made
length(unique(PO_100extra$ID2)) # 95 test f1's assigned
table(table(PO_100extra$ID2)) # 95 assigned to two
length(unique(PO_100extra$ID1)) # 65 of 65 parents assigned
summary(PO_100extra$PO) # very high confidence

PO_100extra_counts <- PO_100extra %>% 
  group_by(ID2) %>%
  mutate(n_parents = n()) %>%
  ungroup()

table(PO_100extra_counts$n_parents) 
PO_100extra_1 <- PO_100extra_counts %>% filter(n_parents == 1) # 0
PO_100extra_2 <- PO_100extra_counts %>% filter(n_parents == 2) # 95
##### ---- ---- #####

### ---- validate the inferred parent/grandparent crosses ---- ###
##### ---- set up cross lookup ---- #####
test_crosses # remember, here's where we have stored our 61 crosses
for(i in 1:nrow(test_crosses)){
  
  p1 <- test_crosses$male[i]
  p2 <- test_crosses$female[i]
  
  id1 <- as.numeric(sub("f0_", "", p1))
  id2 <- as.numeric(sub("f0_", "", p2))
  
  test_crosses$Pair[i] <- paste(
    paste0("f0_", pmin(id1, id2)),
    paste0("f0_", pmax(id1, id2)),
    sep="__"
  )
}

##### ---- best case ---- #####
PO_bestcase_2_valid <- PO_bestcase_2 %>%
  group_by(ID2) %>%
  mutate(
    inferred_pair = {
      parent_nums <- sort(as.numeric(sub("f0_", "", ID1)))
      paste0("f0_", parent_nums, collapse = "__")
    },
    valid_cross = inferred_pair %in% test_crosses$Pair
  ) %>%
  ungroup()
table(PO_bestcase_2_valid$valid_cross)
##### ---- 75% ptps ---- #####
PO_75ptps_2_valid <- PO_75ptps_2 %>%
  group_by(ID2) %>%
  mutate(
    inferred_pair = {
      parent_nums <- sort(as.numeric(sub("f0_", "", ID1)))
      paste0("f0_", parent_nums, collapse = "__")
    },
    valid_cross = inferred_pair %in% test_crosses$Pair
  ) %>%
  ungroup()
table(PO_75ptps_2_valid$valid_cross) ## 88 true means 44 inferred correctly
##### ---- 50% ptps ---- #####
PO_50ptps_2_valid <- PO_50ptps_2 %>%
  group_by(ID2) %>%
  mutate(
    inferred_pair = {
      parent_nums <- sort(as.numeric(sub("f0_", "", ID1)))
      paste0("f0_", parent_nums, collapse = "__")
    },
    valid_cross = inferred_pair %in% test_crosses$Pair
  ) %>%
  ungroup()
table(PO_50ptps_2_valid$valid_cross) ## 38 true means 19 inferred correctly
##### ---- 25% ptps ---- #####
PO_25ptps_2_valid <- PO_25ptps_2 %>%
  group_by(ID2) %>%
  mutate(
    inferred_pair = {
      parent_nums <- sort(as.numeric(sub("f0_", "", ID1)))
      paste0("f0_", parent_nums, collapse = "__")
    },
    valid_cross = inferred_pair %in% test_crosses$Pair
  ) %>%
  ungroup()
table(PO_25ptps_2_valid$valid_cross) ## 10 true means 5 inferred correctly
##### ---- 25% extra ---- #####
PO_25extra_2_valid <- PO_25extra_2 %>%
  group_by(ID2) %>%
  mutate(
    inferred_pair = {
      parent_nums <- sort(as.numeric(sub("f0_", "", ID1)))
      paste0("f0_", parent_nums, collapse = "__")
    },
    valid_cross = inferred_pair %in% test_crosses$Pair
  ) %>%
  ungroup()
table(PO_25extra_2_valid$valid_cross) ## 190 true means 95 inferred correctly
##### ---- 50% extra ---- #####
PO_50extra_2_valid <- PO_50extra_2 %>%
  group_by(ID2) %>%
  mutate(
    inferred_pair = {
      parent_nums <- sort(as.numeric(sub("f0_", "", ID1)))
      paste0("f0_", parent_nums, collapse = "__")
    },
    valid_cross = inferred_pair %in% test_crosses$Pair
  ) %>%
  ungroup()
table(PO_50extra_2_valid$valid_cross) ## 190 true means 95 inferred correctly
##### ---- 75% extra ---- #####
PO_75extra_2_valid <- PO_75extra_2 %>%
  group_by(ID2) %>%
  mutate(
    inferred_pair = {
      parent_nums <- sort(as.numeric(sub("f0_", "", ID1)))
      paste0("f0_", parent_nums, collapse = "__")
    },
    valid_cross = inferred_pair %in% test_crosses$Pair
  ) %>%
  ungroup()
table(PO_75extra_2_valid$valid_cross) ## 190 true means 95 inferred correctly
##### ---- 100% extra ---- #####
PO_100extra_2_valid <- PO_100extra_2 %>%
  group_by(ID2) %>%
  mutate(
    inferred_pair = {
      parent_nums <- sort(as.numeric(sub("f0_", "", ID1)))
      paste0("f0_", parent_nums, collapse = "__")
    },
    valid_cross = inferred_pair %in% test_crosses$Pair
  ) %>%
  ungroup()
table(PO_100extra_2_valid$valid_cross) ## 190 true means 95 inferred correctly
##### ---- ---- #####







