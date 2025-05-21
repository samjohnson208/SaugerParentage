### --- sauger_hiphop_second.R --- ###
### created by SPJ on 050525 ###
# purpose: explore differences in hiphip output with varying filters on input vcfs

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

# soft_genotypes MAF 1, MISS 90. Should perform worse than our original MAF 5, MISS 90.
soft_genotypes <- read.table(file = "soft_variants_bial_noindels_q20_maf1_miss90.012", 
                             header = FALSE, sep = "\t", na.strings = c("NA", "-1"))
soft_genotypes <- soft_genotypes[, -1]
soft_ind <- read.table(file = "soft_variants_bial_noindels_q20_maf1_miss90.012.indv", header = FALSE)
rownames(soft_genotypes) <- soft_ind$V1


hard_genotypes <- read.table(file = "hard_variants_bial_noindels_q20_mindep8_maxdep75_maf5_miss95.012", 
                             header = FALSE, sep = "\t", na.strings = c("NA", "-1"))
hard_genotypes <- hard_genotypes[, -1]
hard_ind <- read.table(file = "hard_variants_bial_noindels_q20_mindep8_maxdep75_maf5_miss95.012.indv", header = FALSE)
rownames(hard_genotypes) <- hard_ind$V1
  
# read in parent offspring (F0 and Test F1 sample id's)
testsamp <- individuals$individual # 208 individuals

#create genotype matrix parent-offspring by filtering rownames of gmmat to include only those in testsamp$sample
soft_genotypes <- soft_genotypes[rownames(soft_genotypes) %in% testsamp, , drop = FALSE]
hard_genotypes <- hard_genotypes[rownames(hard_genotypes) %in% testsamp, , drop = FALSE]
# now that the dataframes have been called, there are formatting changes that need to be made.
  
# ensure that this worked, and that individuals can be treated the same in both filtering cases
table(soft_ind == hard_ind) # 1184 TRUE. nice.

# we can now treat all individuals the same. no need for soft_ind or hard_ind

# now that the dataframes have been called, there are a few formatting changes that need to be made.

# first, individuals$type needs to include the sex of the F0 parents, formatted as "adult male/female"

for(i in 1:nrow(individuals)) {
  if (individuals$type[i] == "M") {
    individuals$type[i] <- "adult male"
  } else if (individuals$type[i] == "F") {
    individuals$type[i] <- "adult female"
  } else if (individuals$type[i] == "offspring") {}
}

# second, genotypes as outputted by vcftools --012 is formatted differently than is required for hiphop
# vcftools:
# 0 = homozygous reference
# 1 = heterozygous reference, alternate
# 2 = homozygous alternate
# -1 = missing

# hiphop:
# 0 = homozygous for one allele
# 1 = homozygous for the alternative allele
# 2 = heterozygous at this locus
# NA = missing

# because we read in genotypes with na.strings = c("NA", "-1"), we're good there,
# but now we need to change all values of genotypes to match the hiphop format
# we need: all 1's to be 2's, all 2's to be 1's.

# apply the following function across all columns of the genotype matrix
# if the value is a 1 turn it to a 2, and vice versa
soft_geno_conv <- as.data.frame(apply(soft_genotypes, 2, function(x) {
  ifelse(x == 1, 2, ifelse(x == 2, 1, x))
}))

hard_geno_conv <- as.data.frame(apply(hard_genotypes, 2, function(x) {
  ifelse(x == 1, 2, ifelse(x == 2, 1, x))
}))

# read in the large sar_data file that houses all of the true crosses and filter
# it to only include the testsamp individuals
sar <- read.csv(file = "SAR_Data_092424.csv", header = TRUE)
sar <- sar %>% 
  filter(Sample_ID %in% testsamp)
sar_clean <- sar %>%
  mutate(Crossed_with_clean = ifelse(Crossed_with == "", NA, Crossed_with)) %>%
  mutate(Crossed_with_clean = strsplit(as.character(Crossed_with), ","))

# alright, now let's go ahead and just remove that individual from all of these dataframes
# and see if our error rate goes down at all.
sar_clean <- sar_clean %>% 
  filter(Sample_ID != "SAR_15_6433")
# down to 207 observations. 95 offspring + 112 parents. 

# use indexing for rownames, filter for columns
soft_geno_conv <- soft_geno_conv[rownames(soft_geno_conv) != "SAR_15_6433", ]
hard_geno_conv <- hard_geno_conv[rownames(hard_geno_conv) != "SAR_15_6433", ]

individuals <- individuals %>% 
  filter(individual != "SAR_15_6433")

# just need to remove it from four dataframes: soft_/hard_geno_conv, individuals, and sar_clean. 
# 207 obs -> inspection

# the inspection done here lets you verify that all of your individuals are present
# in the genotype matrix
soft_inspection <- inspect(ind = individuals, gen = soft_geno_conv)
head(soft_inspection)
soft_inspection[which(soft_inspection$sampled == 0), ] # got everyone!

hard_inspection <- inspect(ind = individuals, gen = hard_geno_conv)
head(hard_inspection)
hard_inspection[which(hard_inspection$sampled == 0), ] # got everyone!

# the inspection done here lets you verify that all of your individuals are present
# in the genotype matrix

# this is what generates for you all of the possible relationships at play
# this is to be subsetted downstream
soft_combinations <- hothiphop(ind = individuals, gen = soft_geno_conv)
hard_combinations <- hothiphop(ind = individuals, gen = hard_geno_conv)

# now we will return the top three most likely relationships for each offspring
# according to the hothiphop score
soft_top3 <- topmatch(x = soft_combinations, ranking = "hothiphop.parents")
soft_top3 <- soft_top3 %>% 
    filter(dam.type != "social parent")

hard_top3 <- topmatch(x = hard_combinations, ranking = "hothiphop.parents")
hard_top3 <- hard_top3 %>% 
  filter(dam.type != "social parent")

# since the lookup table has already been made, just read it in.
cross_lookup <- read.csv(file = "true_cross_lookup.csv", header = TRUE)
cross_lookup <- cross_lookup[,2:4] # came in w/ an extra column... weird.

## JOIN THE LOOKUP TABLE AND THE TOP3 TABLES ##

# this is a ridiculously effective set of piped functions here:
# first, we create this new dataframe from top3
soft_top3_checked <- soft_top3 %>%
  # create a new column called pair, where the cross is structured the same way
  # as those in the lookup table. again, we're using pmin and pmax. sorting the
  # pair entries this way ensures that the cross A_B will be treated the same as
  # the cross B_A.
  mutate(pair = paste(pmin(dam, sire), pmax(dam, sire), sep = "_")) %>%
  
  left_join(cross_lookup %>% # join the top3 dataframe to the lookup table
              select(pair) %>% # but ONLY the column cross_lookup$pair
              distinct() %>% # keeps only unique entries of the known pairs, avoids any repeats
              mutate(valid_cross = TRUE), # creats a new column called valid_cross, which is either TRUE
            # or NA, if the pair exists in top3, or if it does not.
            by = "pair") %>% # actually completes the join, by the shared pair column.
  # this is the line in which the matching actually occurs. 
  mutate(valid_cross = ifelse(is.na(valid_cross), FALSE, valid_cross)) # recreates that valid_cross column but turns
# all of the NA's formed by the join into FALSE's

hard_top3_checked <- hard_top3 %>%
  # create a new column called pair, where the cross is structured the same way
  # as those in the lookup table. again, we're using pmin and pmax. sorting the
  # pair entries this way ensures that the cross A_B will be treated the same as
  # the cross B_A.
  mutate(pair = paste(pmin(dam, sire), pmax(dam, sire), sep = "_")) %>%
  
  left_join(cross_lookup %>% # join the top3 dataframe to the lookup table
              select(pair) %>% # but ONLY the column cross_lookup$pair
              distinct() %>% # keeps only unique entries of the known pairs, avoids any repeats
              mutate(valid_cross = TRUE), # creats a new column called valid_cross, which is either TRUE
            # or NA, if the pair exists in top3, or if it does not.
            by = "pair") %>% # actually completes the join, by the shared pair column.
  # this is the line in which the matching actually occurs. 
  mutate(valid_cross = ifelse(is.na(valid_cross), FALSE, valid_cross)) # recreates that valid_cross column but turns
# all of the NA's formed by the join into FALSE's

# after the join, the columns of interest are all on the right side of the dataframe, out of sight.
# rearrange using select
hard_top3_checked <- hard_top3_checked %>%
  select(year, brood, offspring, rank, dam, sire, pair, valid_cross, everything())

soft_top3_checked <- soft_top3_checked %>% 
  select(year, brood, offspring, rank, dam, sire, pair, valid_cross, everything())  

##### ----- Some summary information ----- #####
## -- Soft Filtering -- ##
length(unique(soft_top3_checked$offspring)) # for 95 test F1's
table(soft_top3_checked$valid_cross) # we're capturing 88 true relationships by using 
# the top 3 ranked pairs for each offspring

# among the top ranked parent pair for each offspring...
soft_rank1 <- soft_top3_checked %>% 
  filter(rank == 1)
table(soft_rank1$valid_cross) # we captured the correct relationship 85/95 times 
# (88.42% accuracy)

# among the second highest ranked parent pair for each offspring...
soft_rank2 <- soft_top3_checked %>% 
  filter(rank == 2)
table(soft_rank2$valid_cross) # we captured the correct relationship 0 times...

# and among the third highest ranked parent pair for each offspring...
soft_rank3 <- soft_top3_checked %>% 
  filter(rank == 3) 
table(soft_rank3$valid_cross) # we captured the correct relationship an additional 3 times..

# WEIRD THAT WE LOSE ALL OF THE SECOND MOST LIKELY RELATIONSHIPS 


## -- Hard Filtering -- ##
length(unique(hard_top3_checked$offspring)) # for 95 test F1's
table(hard_top3_checked$valid_cross) # we're capturing 92 true relationships by using 
# the top 3 ranked pairs for each offspring

# among the top ranked parent pair for each offspring...
hard_rank1 <- hard_top3_checked %>% 
  filter(rank == 1)
table(hard_rank1$valid_cross) # we captured the correct relationship 85/95 times 
# (88.42% accuracy)

# among the second highest ranked parent pair for each offspring...
hard_rank2 <- hard_top3_checked %>% 
  filter(rank == 2)
table(hard_rank2$valid_cross) # we captured the correct relationship an additional 6 times...

# and among the third highest ranked parent pair for each offspring...
hard_rank3 <- hard_top3_checked %>% 
  filter(rank == 3) 
table(hard_rank3$valid_cross) # we captured the correct relationship an additional 1 time...

# Now I want to see if the percentage of missing data per individual is indicative of our
# ability to generate the correct relationships.

# do the rownames of each genotype matrix agree?
table(rownames(soft_geno_conv) == rownames(hard_geno_conv))

md_per_individual <- data.frame(rownames(soft_geno_conv), soft_md = NA, hard_md = NA)
md_per_individual$soft_md <- rowMeans(is.na(soft_geno_conv)) * 100
md_per_individual$hard_md <- rowMeans(is.na(hard_geno_conv)) * 100
summary(md_per_individual$soft_md)
summary(md_per_individual$hard_md)

# Alright. Here's what I'm gathering. 
# Focusing primarily on In the cases where either:
  # 1. the top ranked relationship is false
  # 2. the second ranked relationship is true
  # 3. the third ranked relationship is true

# It's the same few parents that are being picked out by hiphop. 

# Message to Katie and Josh:
# Stricter filtering gets us to 92/95 true parent pairs returned! 
# (85/95 from the top ranked/most likely relationship for each offspring). 
# MD per individual does not appear to be driving whether or not the true relationship 
# is returned. Seems like in cases where the top ranked relationship is false, 
# or the second or third ranked relationship is true, it’s the same few parents 
# that are being picked out.

# An example: One of the true crosses is SAR_15_6401 x SAR_15_6417. That’s returned 
# as true in a few of the second ranked relationships, but not in the top ranked 
# relationships. In the top ranked relationships, SAR_15_6401 x SAR_15_6440 often appears. 
# Maybe SAR_15_6417 and SAR_15_6440 are siblings or closely related in some other way?

# Filtering for this hiphop run was for only biallelic sites, no indels, phred > 20, 
# mindepth > 8, maxdepth < 75, MAF >= 5%, Missing Data Per Site =< 5%. (HARD FILTER)

# Number of sites = 4196

# I’m still unsure as to how my number of sites after this new filtering approach 
# was so much higher than what I was putting into Sequoia… I’ll have to investigate this…
# That sequoia dataset was MAF 30, Miss 90, that's why.

# Let's inspect the hothiphop.parents scores for the true and false relationships
hard_rank1false <- hard_rank1 %>% 
  filter(valid_cross == FALSE) # want the false top ranked relationships
hard_rank2true <- hard_rank2 %>% 
  filter(valid_cross == TRUE) # and the true second...
hard_rank3true <- hard_rank3 %>% 
  filter(valid_cross == TRUE) # and third relationships...


tfpositions <- as.numeric(factor(hard_top3_checked$valid_cross))
tfpositionshardrank2true <- c(2, 2, 2, 2, 2, 2)
tfpositionshardrank3true <- 2
tfpositionshardrank1false <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)


boxplot(hard_top3_checked$hothiphop.parents ~ hard_top3_checked$valid_cross, outline = FALSE)

# the goal is to remove duplicate points on these plots. as of now, we'll double plot instances
# where the top ranked relationship is false, and instances where the second and third are true.
# this compound filter should remove those.
hard_to_plot <- hard_top3_checked %>% 
    filter(!(rank == 1 & valid_cross == FALSE)) %>% 
    filter(!(rank == 2 & valid_cross == TRUE)) %>% 
    filter(!(rank == 3 & valid_cross == TRUE))

hard_to_plot <- data.frame(hard_to_plot, tfpositions = NA)
  
# tfpositionshardtoplot <- as.numeric(factor(hard_to_plot$valid_cross))
hard_to_plot$tfpositions <- as.numeric(factor(hard_to_plot$valid_cross))
  

boxplot(hard_top3_checked$hothiphop.parents ~ hard_top3_checked$valid_cross, outline = FALSE, 
        main = "F0/Test F1 Relationships", xlab = "Validation of Inferred Cross", 
        ylab = "Parent Pair HotHipHop Score")
points(jitter(tfpositionshardtoplot, amount = 0.125), hard_to_plot$hothiphop.parents, pch = 19, col = tfpositionshardtoplot + 1)
# points(jitter(tfpositions, amount = 0.15), hard_top3_checked$hothiphop.parents, pch = 19, col = tfpositions + 1)
points(jitter(tfpositionshardrank2true, amount = 0.02), hard_rank2true$hothiphop.parents, pch = 19, col = "darkgreen")
points(jitter(tfpositionshardrank3true, amount = 0.015), hard_rank3true$hothiphop.parents, pch = 19, col = "darkgreen")
points(jitter(tfpositionshardrank1false, amount = 0.02), hard_rank1false$hothiphop.parents, pch = 19, col = "darkred")
legend("bottomleft", legend = c("1st Rank - TRUE", "2nd/3rd Rank - TRUE", "1st Rank - FALSE", "2nd/3rd Rank - FALSE"),
       col = c(3, "darkgreen", "darkred", 2), pch = 19, cex = 1.1)
abline(h = 230, lty = 2, lwd = 2)

write.csv(hard_geno_conv, "hardfilter_genotypes_conv.csv")
write.csv(hard_top3_checked, "hardfilter_top3_checked.csv")
write.csv(hard_to_plot, "hothiphop_hard_to_plot_100.csv")

