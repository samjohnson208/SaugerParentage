### --- sauger_hiphop_third.R --- ###
### created by SPJ on 050725 ###
# purpose: explore differences in hiphip output with varying PTPS for F0's
# test run for 50% ptps

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
      set.seed(12)
      hard_geno_conv_F0_50 <- hard_geno_conv_F0[sample(nrow(hard_geno_conv_F0), size = 56),]
      # random 56 (50%) of the parents!

      # filter for f1's
      hard_geno_conv_test_f1 <- hard_geno_conv %>% 
        filter(grepl("^SAR_15_67", rownames(hard_geno_conv)))

# recombine to create genotype matrix upon which hiphop will act
hard_geno_conv_50 <- rbind(hard_geno_conv_F0_50, hard_geno_conv_test_f1)

# as well as individuals dataframe for this subset of f0's and all f1's
individuals_50 <- individuals %>% 
    filter(individual %in% rownames(hard_geno_conv_50))

# prepare indivs using hiphop's specs
for(i in 1:nrow(individuals_50)) {
  if (individuals_50$type[i] == "M") {
    individuals_50$type[i] <- "adult male"
  } else if (individuals_50$type[i] == "F") {
    individuals_50$type[i] <- "adult female"
  } else if (individuals_50$type[i] == "offspring") {}
}

# make sure we have everyone we want to account for
inspection_50 <- inspect(ind = individuals_50, gen = hard_geno_conv_50)
inspection_50[which(inspection_50$sampled == 0), ] #consistent. nice.

# run hothiphop
combinations50 <- hothiphop(ind = individuals_50, gen = hard_geno_conv_50)

# keep the top 3 parent relationships for each offspring, remove social parent rows
top3_50 <- topmatch(x = combinations50, ranking = "hothiphop.parents")
top3_50 <- top3_50 %>% 
    filter(dam.type != "social parent")

# filter for only top ranked relationships
top3_50_rank1 <- top3_50 %>%
    filter(rank == 1)

# now generates the inferred pair column for this ptps
top3_50_rank1 <- top3_50_rank1 %>% 
  mutate(pair = paste(pmin(dam, sire), pmax(dam, sire), sep = "_"))

# prepare this dataframe to merge with the 100 ptps observations
top3_50_rank1_tomerge <- data.frame(offspring = top3_50_rank1$offspring,
                                    dam_50 = top3_50_rank1$dam, 
                                    sire_50 = top3_50_rank1$sire,
                                    pair_50 = top3_50_rank1$pair,
                                    score_50 = top3_50_rank1$hothiphop.parents,
                                    valid_cross_50 = NA)





# read in the observed relationships for 100% ptps
hard_top3_checked <- read.csv(file = "hardfilter_top3_checked.csv", header = TRUE)
hard_top3_checked <- hard_top3_checked[, -1]
    # keep only the top ranked relationships, and only those that are determined to be true (85/95)
    # for the ptps analysis to work, we need to assume that these are all true relationships.
    # these relationships will now serve as the lookup table for additional ptps situations.
      hard_top3_checked_rank1 <- hard_top3_checked %>% 
        filter(rank == 1)
      hard_top3_checked_rank1_true <- hard_top3_checked_rank1 %>% 
        filter(valid_cross == TRUE)

# Now 85 offspring with true top ranked relationships.
# So now we need to start making our new dataframe.

hard_top3_100 <- data.frame(brood = hard_top3_checked_rank1_true$brood, 
                            offspring = hard_top3_checked_rank1_true$offspring, 
                            dam_100 = hard_top3_checked_rank1_true$dam, 
                            sire_100 = hard_top3_checked_rank1_true$sire,
                            pair_100 = hard_top3_checked_rank1_true$pair,
                            valid_cross = hard_top3_checked_rank1_true$valid_cross, 
                            score_100 = hard_top3_checked_rank1_true$hothiphop.parents)

# again, make sure that the offspring that are going into the merged df are the same across both ptps situations
top3_50_rank1_tomerge <- top3_50_rank1_tomerge %>% 
    filter(top3_50_rank1_tomerge$offspring %in% hard_top3_100$offspring)

# merge the dfs. same information now for both ptps situations. 
top1_100_50 <- merge(hard_top3_100, top3_50_rank1_tomerge, by = "offspring")

# now see if the pair given 50% ptps is the same as 100% ptps, and validate it in the valid_cross_50 column
for(i in 1:nrow(top1_100_50)){
  p_50 <- top1_100_50$pair_50[i]
  p_100 <- top1_100_50$pair_100[i]
  top1_100_50$valid_cross_50[i] <- ifelse(p_50 == p_100, TRUE, FALSE)
}


table(top1_100_50$valid_cross_50)
# With 50% PTPS, we recover the true relationship 41% of the time.


# set up a count column to see what the hothiphop scores look like for 50% ptps when each offspring has 0, 1, 2 parents sampled.
top1_100_50 <- data.frame(top1_100_50, count_par_50 = 0)

# count how many of the true parents for an individual were included in this 50% subset
for(i in 1:nrow(top1_100_50)) {
  if(top1_100_50$dam_100[i] %in% rownames(hard_geno_conv_F0_50)) {
    top1_100_50$count_par_50[i] <- top1_100_50$count_par_50[i] + 1
  }
  if(top1_100_50$sire_100[i] %in% rownames(hard_geno_conv_F0_50)) {
    top1_100_50$count_par_50[i] <- top1_100_50$count_par_50[i] + 1
  }
}
table(top1_100_50$count_par_50)
# With 50% PTPS, we recover the true relationship 41% of the time.
# However, when we don't get both parents, we pick up 1 true parent another 51% of the time.
# Therefore, there is only 9% of the time that we don't get either true parent.

# set up positions for plotting
top1_100_50$pspositions <- as.numeric(factor(top1_100_50$count_par_50 + 1))

# filter for adding points
ps50_0 <- top1_100_50 %>% 
    filter(count_par_50 == 0)
ps50_0$pspositions <- 1
ps50_1 <- top1_100_50 %>% 
  filter(count_par_50 == 1)
ps50_1$pspositions <- 2
ps50_2 <- top1_100_50 %>% 
  filter(count_par_50 == 2)
ps50_2$pspositions <- 3


boxplot(top1_100_50$score_50 ~ top1_100_50$count_par_50, outline = FALSE, 
        main = "50% True Parents Sampled", xlab = "Number of Parents Sampled Per Offspring", 
        ylab = "Parent Pair HotHipHop Score")
points(jitter(ps50_0$pspositions, amount = 0.125), ps50_0$score_50, pch = 19, col = "darkred")
points(jitter(ps50_1$pspositions, amount = 0.125), ps50_1$score_50, pch = 19, col = "darkred")
points(jitter(ps50_2$pspositions, amount = 0.125), ps50_2$score_50, pch = 19, col = "darkgreen")

write_csv(top1_100_50, "top1_100_50.csv")







