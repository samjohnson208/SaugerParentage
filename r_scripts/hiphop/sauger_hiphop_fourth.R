### --- sauger_hiphop_fourth.R --- ###
### created by SPJ on 050725 ###
# purpose: explore differences in hiphip output with varying PTPS for F0's
# add 25% and 75% PTPS

getwd()
setwd("/Users/samjohnson/Desktop/hiphop")

#install.packages("hiphop")
library(hiphop)
library(dplyr)
library(tidyverse)
library(purrr)

# read in a an individuals dataframe as formatted according to the hiphop documentation
individuals <- read.csv(file = "testindivs.csv", header = TRUE)
# 95 test f1's, 113 F0 parents

# read in genotypes
hard_geno_conv <- read.csv("hardfilter_genotypes_conv.csv", header = TRUE)
rownames(hard_geno_conv) <- hard_geno_conv$X
hard_geno_conv <- hard_geno_conv[, -1]

# separate into f0's and f1's
hard_geno_conv_F0 <- hard_geno_conv %>% 
  filter(grepl("^SAR_15_64|^SAR_15_65", rownames(hard_geno_conv)))
# 112 Parents. All good.

      # take a random sample of the f0's
      set.seed(720)
      hard_geno_conv_F0_75 <- hard_geno_conv_F0[sample(nrow(hard_geno_conv_F0), size = 84),]
      # random 84 (75%) of the parents!
      
      # filter for f1's
      hard_geno_conv_test_f1 <- hard_geno_conv %>% 
        filter(grepl("^SAR_15_67", rownames(hard_geno_conv))) #95 f1's

# recombine to create genotype matrix upon which hiphop will act
hard_geno_conv_75 <- rbind(hard_geno_conv_F0_75, hard_geno_conv_test_f1) # 179 individuals

# as well as individuals dataframe for this subset of f0's and all f1's
individuals_75 <- individuals %>% 
  filter(individual %in% rownames(hard_geno_conv_75))

# prepare indivs using hiphop's specs
for(i in 1:nrow(individuals_75)) {
  if (individuals_75$type[i] == "M") {
    individuals_75$type[i] <- "adult male"
  } else if (individuals_75$type[i] == "F") {
    individuals_75$type[i] <- "adult female"
  } else if (individuals_75$type[i] == "offspring") {}
}

# make sure we have everyone we want to account for
inspection_75 <- inspect(ind = individuals_75, gen = hard_geno_conv_75)
inspection_75[which(inspection_75$sampled == 0), ] #consistent. nice.

# run hothiphop
combinations75 <- hothiphop(ind = individuals_75, gen = hard_geno_conv_75)

# keep the top 3 parent relationships for each offspring, remove social parent rows
top3_75 <- topmatch(x = combinations75, ranking = "hothiphop.parents")
top3_75 <- top3_75 %>% 
  filter(dam.type != "social parent")

# filter for only top ranked relationships
top3_75_rank1 <- top3_75 %>%
  filter(rank == 1)
      
# now generates the inferred pair column for this ptps
top3_75_rank1 <- top3_75_rank1 %>% 
  mutate(pair = paste(pmin(dam, sire), pmax(dam, sire), sep = "_"))
      
# prepare this dataframe to merge with the 100 ptps observations
top3_75_rank1_tomerge <- data.frame(offspring = top3_75_rank1$offspring,
                                    dam_75 = top3_75_rank1$dam, 
                                    sire_75 = top3_75_rank1$sire,
                                    pair_75 = top3_75_rank1$pair,
                                    score_75 = top3_75_rank1$hothiphop.parents,
                                    valid_cross_75 = NA)   

# read in the observed relationships for 100 and 50% ptps      
top1_100_50 <- read.csv(file = "top1_100_50.csv", header = TRUE)

# just be safe and filter using this method instead of letting merge delete rows
top3_75_rank1_tomerge <- top3_75_rank1_tomerge %>% 
  filter(top3_75_rank1_tomerge$offspring %in% top1_100_50$offspring)

# merge the dfs. same information now for both ptps situations. 
top1_100_75_50 <- merge(top1_100_50, top3_75_rank1_tomerge, by = "offspring")

# now see if the pair given 75% ptps is the same as 100% ptps, and validate it in the valid_cross_75 column
for(i in 1:nrow(top1_100_75_50)){
  p_75 <- top1_100_75_50$pair_75[i]
  p_100 <- top1_100_75_50$pair_100[i]
  top1_100_75_50$valid_cross_75[i] <- ifelse(p_75 == p_100, TRUE, FALSE)
}

table(top1_100_75_50$valid_cross_75)
# With 75% PTPS, we recover the true relationship 74% of the time.

# set up a count column to see what the hothiphop scores look like for 75% ptps when each offspring has 0, 1, 2 parents sampled.
top1_100_75_50 <- data.frame(top1_100_75_50, count_par_75 = 0)

# count how many of the true parents for an individual were included in this 50% subset
for(i in 1:nrow(top1_100_75_50)) {
  if(top1_100_75_50$dam_100[i] %in% rownames(hard_geno_conv_F0_75)) {
    top1_100_75_50$count_par_75[i] <- top1_100_75_50$count_par_75[i] + 1
  }
  if(top1_100_75_50$sire_100[i] %in% rownames(hard_geno_conv_F0_75)) {
    top1_100_75_50$count_par_75[i] <- top1_100_75_50$count_par_75[i] + 1
  }
}
table(top1_100_75_50$count_par_75)
# With 75% PTPS, we recover the true relationship 43% of the time.
# However, when we don't get both parents, we pick up 1 true parent another 38% of the time.
# Therefore, we only assign two incorrect parents 17% of the time, heh. (for THIS REPLICATE)

# set up positions for plotting. had to rename the ps positions from the 50% ptps run.
top1_100_75_50 <- top1_100_75_50 %>% 
    rename(pspositions_50 = pspositions)

top1_100_75_50$pspositions_75 <- as.numeric(factor(top1_100_75_50$count_par_75 + 1))

# filter for adding points
ps75_0 <- top1_100_75_50 %>% 
  filter(count_par_75 == 0)
ps75_1 <- top1_100_75_50 %>% 
  filter(count_par_75 == 1)
ps75_2 <- top1_100_75_50 %>% 
  filter(count_par_75 == 2)



boxplot(top1_100_75_50$score_75 ~ top1_100_75_50$count_par_75, outline = FALSE, 
        main = "75% True Parents Sampled", xlab = "Number of Parents Sampled Per Offspring", 
        ylab = "Parent Pair HotHipHop Score")
points(jitter(ps75_0$pspositions_75, amount = 0.125), ps75_0$score_75, pch = 19, col = "darkred")
points(jitter(ps75_1$pspositions_75, amount = 0.125), ps75_1$score_75, pch = 19, col = "darkred")
points(jitter(ps75_2$pspositions_75, amount = 0.125), ps75_2$score_75, pch = 19, col = "darkgreen")

write_csv(top1_100_75_50, "top1_100_75_50.csv")

# Alright. This is absolutely stupid, and here is why. 
# These results don't mean ANYTHING without replication. If you change the seed
# and pull a DIFFERENT 75% of the parents, your results change. For example, I 
# also tried this with a set.seed(12), and I had NO instances where I had 0
# parents sampled for any of the offspring. Think this is going to have to be
# replicated LOTS of times to see where these distributions truly fall. This is
# the right idea, but it doesn't hold any water without replication.

# I'll come back at some point and add all of this code for 25% ptps also, but
# we're going to have to come up with an absolutely mind boggling script that will
# iterate this WHOLE thing (sampling parents -> running hiphop -> summarizing dist of scores)
# a bunch of times... obviously not going to happen for this talk. That's a week long project
# just to get the script together, let alone to run it, analyze, and generate figures...

setwd("/Users/samjohnson/Desktop/hiphop/")
top1_100_75_50 <- read.csv(file = "top1_100_75_50.csv", header = TRUE)

colnames(top1_100_75_50)

table(top1_100_75_50$sire_100)
table(top1_100_75_50$dam_100)

table(top1_100_75_50$sire_75)
table(top1_100_75_50$dam_75)

table(top1_100_75_50$sire_50)
table(top1_100_75_50$dam_50)








#### ---- CHECK FOR 25% PARENTS SAMPLED ---- ####

# read in a an individuals dataframe as formatted according to the hiphop documentation
individuals <- read.csv(file = "testindivs.csv", header = TRUE)
# 95 test f1's, 113 F0 parents

# read in genotypes
hard_geno_conv <- read.csv("hardfilter_genotypes_conv.csv", header = TRUE)
rownames(hard_geno_conv) <- hard_geno_conv$X
hard_geno_conv <- hard_geno_conv[, -1]

# separate into f0's and f1's
hard_geno_conv_F0 <- hard_geno_conv %>% 
  filter(grepl("^SAR_15_64|^SAR_15_65", rownames(hard_geno_conv)))
# 112 Parents. All good.

# take a random sample of the f0's
set.seed(720)
hard_geno_conv_F0_25 <- hard_geno_conv_F0[sample(nrow(hard_geno_conv_F0), size = 28),]
# random 84 (75%) of the parents!

# filter for f1's
hard_geno_conv_test_f1 <- hard_geno_conv %>% 
  filter(grepl("^SAR_15_67", rownames(hard_geno_conv))) #95 f1's

# recombine to create genotype matrix upon which hiphop will act
hard_geno_conv_25 <- rbind(hard_geno_conv_F0_25, hard_geno_conv_test_f1) # 123 individuals

# as well as individuals dataframe for this subset of f0's and all f1's
individuals_25 <- individuals %>% 
  filter(individual %in% rownames(hard_geno_conv_25))

# prepare indivs using hiphop's specs
for(i in 1:nrow(individuals_25)) {
  if (individuals_25$type[i] == "M") {
    individuals_25$type[i] <- "adult male"
  } else if (individuals_25$type[i] == "F") {
    individuals_25$type[i] <- "adult female"
  } else if (individuals_25$type[i] == "offspring") {}
}

# make sure we have everyone we want to account for
inspection_25 <- inspect(ind = individuals_25, gen = hard_geno_conv_25)
inspection_25[which(inspection_25$sampled == 0), ] #consistent. nice.

# run hothiphop
combinations25 <- hothiphop(ind = individuals_25, gen = hard_geno_conv_25)

# keep the top 3 parent relationships for each offspring, remove social parent rows
top3_25 <- topmatch(x = combinations25, ranking = "hothiphop.parents")
top3_25 <- top3_25 %>% 
  filter(dam.type != "social parent")

# filter for only top ranked relationships
top3_25_rank1 <- top3_25 %>%
  filter(rank == 1)

# now generates the inferred pair column for this ptps
top3_25_rank1 <- top3_25_rank1 %>% 
  mutate(pair = paste(pmin(dam, sire), pmax(dam, sire), sep = "_"))

# prepare this dataframe to merge with the 100 ptps observations
top3_25_rank1_tomerge <- data.frame(offspring = top3_25_rank1$offspring,
                                    dam_25 = top3_25_rank1$dam, 
                                    sire_25 = top3_25_rank1$sire,
                                    pair_25 = top3_25_rank1$pair,
                                    score_25 = top3_25_rank1$hothiphop.parents,
                                    valid_cross_25 = NA)   

# read in the observed relationships for 100, 75, and 50% ptps      
top1_100_75_50 <- read.csv(file = "top1_100_75_50.csv", header = TRUE)

# just be safe and filter using this method instead of letting merge delete rows
top3_25_rank1_tomerge <- top3_25_rank1_tomerge %>% 
  filter(top3_25_rank1_tomerge$offspring %in% top1_100_75_50$offspring)

# merge the dfs. same information now for both ptps situations. 
top1_100_75_50_25 <- merge(top1_100_50, top3_25_rank1_tomerge, by = "offspring")


# now see if the pair given 25% ptps is the same as 100% ptps, and validate it in the valid_cross_25 column
for(i in 1:nrow(top1_100_75_50_25)){
  p_25 <- top1_100_75_50_25$pair_25[i]
  p_100 <- top1_100_75_50_25$pair_100[i]
  top1_100_75_50_25$valid_cross_25[i] <- ifelse(p_25 == p_100, TRUE, FALSE)
}

table(top1_100_75_50_25$valid_cross_25)
# With 75% PTPS, we recover the true relationship 74% of the time.

# set up a count column to see what the hothiphop scores look like for 25% ptps when each offspring has 0, 1, 2 parents sampled.
top1_100_75_50_25 <- data.frame(top1_100_75_50_25, count_par_25 = 0)

# count how many of the true parents for an individual were included in this 50% subset
for(i in 1:nrow(top1_100_75_50_25)) {
  if(top1_100_75_50_25$dam_100[i] %in% rownames(hard_geno_conv_F0_25)) {
    top1_100_75_50_25$count_par_25[i] <- top1_100_75_50_25$count_par_25[i] + 1
  }
  if(top1_100_75_50_25$sire_100[i] %in% rownames(hard_geno_conv_F0_25)) {
    top1_100_75_50_25$count_par_25[i] <- top1_100_75_50_25$count_par_25[i] + 1
  }
}
table(top1_100_75_50_25$count_par_25)
# With 25% PTPS, we recover the true relationship 5% of the time.
# However, when we don't get both parents, we pick up 1 true parent another 22% of the time.
# Therefore, we only assign two incorrect parents 73% of the time, heh. (for THIS REPLICATE)

# set up positions for plotting.
top1_100_75_50_25$pspositions_25 <- as.numeric(factor(top1_100_75_50_25$count_par_25 + 1))

# filter for adding points
ps25_0 <- top1_100_75_50_25 %>% 
  filter(count_par_25 == 0)
ps25_1 <- top1_100_75_50_25 %>% 
  filter(count_par_25 == 1)
ps25_2 <- top1_100_75_50_25 %>% 
  filter(count_par_25 == 2)

boxplot(top1_100_75_50_25$score_25 ~ top1_100_75_50_25$count_par_25, outline = FALSE, 
        main = "25% True Parents Sampled", xlab = "Number of Parents Sampled Per Offspring", 
        ylab = "Parent Pair HotHipHop Score")
points(jitter(ps25_0$pspositions_25, amount = 0.125), ps25_0$score_25, pch = 19, col = "darkred")
points(jitter(ps25_1$pspositions_25, amount = 0.125), ps25_1$score_25, pch = 19, col = "darkred")
points(jitter(ps25_2$pspositions_25, amount = 0.125), ps25_2$score_25, pch = 19, col = "darkgreen")

write_csv(top1_100_75_50_25, "top1_100_75_50_25.csv")

table(top1_100_75_50_25$sire_100)
table(top1_100_75_50_25$dam_100)

table(top1_100_75_50_25$sire_75)
table(top1_100_75_50_25$dam_75)

table(top1_100_75_50_25$sire_50)
table(top1_100_75_50_25$dam_50)

table(top1_100_75_50_25$sire_25)
table(top1_100_75_50_25$dam_25)


# Lots of instances where we have 0 of the true parents sampled. distributions still look alright though...







