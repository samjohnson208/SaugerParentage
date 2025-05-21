### --- sauger_hiphop_first.R --- ###
### created by SPJ on 042525 ###
# purpose: explore hiphop parentage assignment method and run it on f0's + test-f1's 

getwd()
setwd("/Users/samjohnson/Desktop/")

install.packages("hiphop")
library(hiphop)
library(dplyr)
library(tidyverse)
library(purrr)

# read in a an individuals dataframe as formatted according to the hiphop documentation
individuals <- read.csv(file = "testindivs.csv", header = TRUE)
# 95 test f1's, 113 F0 parents

# read in the genotype matrix from the filtered data
genotypes <- read.table(file = "variants_maf5_miss9.012", header = FALSE, sep = "\t", 
                  na.strings = c("NA", "-1"))

# correct row names as sample id's
genotypes <- genotypes[, -1]
ind <- read.table(file = "variants_maf5_miss9.012.indv", header = FALSE)
rownames(genotypes) <- ind$V1

# read in parent offspring (F0 and Test F1 sample id's)
testsamp <- individuals$individual # 208 individuals

#create genotype matrix parent-offspring by filtering rownames of gmmat to include only those in testsamp$sample
genotypes <- genotypes[rownames(genotypes) %in% testsamp, , drop = FALSE]

# now that the dataframe have been called, there are formatting changes that need to be made.

# first, individuals$type needs to include the sex of the F0 parents
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

        # 0 = homozygous for one allele
        # 1 = homozygous for the alternative allele
        # 2 = heterozygous at this locus
        # NA = missing

# because we read in genotypes with na.strings = c("NA", "-1"), we're good there,
# but now we need to change all values of genotypes to match the hiphop format
# we need: all 1's to be 2's, all 2's to be 1's.

# apply the following function across all columns of the genotype matrix
# if the value is a 1 turn it to a 2, and vice versa
geno_conv <- as.data.frame(apply(genotypes, 2, function(x) {
  ifelse(x == 1, 2, ifelse(x == 2, 1, x))
}))

# the inspection done here lets you verify that all of your individuals are present
# in the genotype matrix
inspection <- inspect(ind = individuals, gen = geno_conv)
head(inspection)
inspection[which(inspection$sampled == 0), ] # got everyone!

# this is what generates for you all of the possible relationships at play
# this is to be subsetted downstream
combinations <- hothiphop(ind = individuals, gen = geno_conv)

# here is a histogram of the hothiphopscores for all relationships
# not entirely informative, just to look
par(mfrow = c(1,1))
hist(combinations$hothiphop.parents)

# now we will return the top three most likely relationships for each offspring
# according to the hothiphop score
top3 <- topmatch(x = combinations, ranking = "hothiphop.parents")







## -- Exploration of the inferred relationships -- ##
# remove social parent entries (all NA anyway in this case)
top3 <- top3 %>% 
    filter(dam.type != "social parent")

# read in the large sar_data file that houses all of the true crosses and filter
# it to only include the testsamp individuals
sar <- read.csv(file = "SAR_Data_092424.csv", header = TRUE)
sar <- sar %>% 
filter(Sample_ID %in% testsamp)

# this neat chunk of code does two important things. 
# 1. if, for an individual, the Crossed_with entry is empty, turn it to an NA
# 2. take each Crossed_with entry, and if it consists of two individuals, split them
#    by the comma that separates the two sample id's. this prepares you for the
#    creation of the lookup table below
sar_clean <- sar %>%
  mutate(Crossed_with_clean = ifelse(Crossed_with == "", NA, Crossed_with)) %>%
  mutate(Crossed_with_clean = strsplit(as.character(Crossed_with), ","))

# now we need to create a lookup table against which we will check all of the 
# mother-father relationships inferred by topmatch() and stored in top3
cross_lookup <- sar_clean %>%
  # keep only each sample, and who it's crossed with
  select(Sample_ID, Crossed_with_clean) %>%
  # this is fascinating here. this takes each of the instances where there are 3 
  # or more individuals included in a cross, and explodes them into their own row
  # of this lookup table. now instead of having a single row describe a cross between
  # 3 individuals, it expands to two (or three) rows of pairwise crosses
  # it ALSO creates a pair entry for each cross in the lookup table, where the cross is listed
  # as being sorted where the smaller sample id is first, so it solves any ordering issues
  unnest(Crossed_with_clean) %>%
  rename(partner = Crossed_with_clean) %>%
  filter(!is.na(partner)) %>%
  # this last line creates the pairwise cross, sorted so that the alphabetically 
  # smaller sample id appears first, and the two id's are separated by an underscore
  mutate(pair = paste(pmin(Sample_ID, partner), pmax(Sample_ID, partner), sep = "_"))

# some of the original entries in the sar_data dataframe included spaces.
# this code tidies that up in the lookup table using str_replace_all(spaces with no spaces)
# and does so in all columns.
library(stringr)
cross_lookup <- cross_lookup %>%
  mutate(
    Sample_ID = str_replace_all(str_trim(Sample_ID), " ", ""),
    partner = str_replace_all(str_trim(partner), " ", ""),
    pair = paste(pmin(Sample_ID, partner), pmax(Sample_ID, partner), sep = "_")
  )

# this is a ridiculously effective set of piped functions here:
# first, we create this new dataframe from top3
top3_checked <- top3 %>%
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
top3_checked <- top3_checked %>%
  select(year, brood, offspring, rank, dam, sire, pair, valid_cross, everything())

## -- Some summary information -- ##
length(unique(top3_checked$offspring)) # for 95 test F1's
table(top3_checked$valid_cross) # we're capturing 90 true relationships by using 
                                # the top 3 ranked pairs for each offspring

# among the top ranked parent pair for each offspring...
rank1 <- top3_checked %>% 
    filter(rank == 1)
table(rank1$valid_cross) # we captured the correct relationship 84/95 times 
                         # (88.42% accuracy)

# among the second highest ranked parent pair for each offspring...
rank2 <- top3_checked %>% 
    filter(rank == 2)
table(rank2$valid_cross) # we captured the correct relationship an additional 3 times...

# and among the third highest ranked parent pair for each offspring...
rank3 <- top3_checked %>% 
  filter(rank == 3)
table(rank3$valid_cross) # we captured the correct relationship an additional 3 times..

# So for our 90/95 correct relationships...
# 84 came from the highest ranked pairs,
# 3 came from the second highest ranked pairs, and
# 3 came from the third highest ranked pairs.

# My question is... what happened to the other 5 offspring? Was it that:
# a) their true parent pairs are ranked lower than the third most likely?
# b) the hothiphop exclusion method didn't have the ability to recognize them at all?
# c) their parents were crossed but not documented?

# we also have a weird case with SAR_15_6433, where it was not documented to have
# been crossed with any other F0's. We should rerun this without including that
# individual. it is (correctly) not in the lookup table, but the whole point here
# is to ONLY give the program the opportunity to assign relationships that we know
# are possible given the cross information that we were given by WGFD.















#### ---- 04/30/25 Work ---- ####
# we have one individual, SAR_15_6433, that was not crossed according to sar_data

getwd()
setwd("/Users/samjohnson/Desktop/hiphop/")

install.packages("hiphop")
library(hiphop)
library(dplyr)
library(tidyverse)
library(purrr)

# read in a an individuals dataframe as formatted according to the hiphop documentation
individuals <- read.csv(file = "testindivs.csv", header = TRUE)
# 95 test f1's, 113 F0 parents

# read in the genotype matrix from the filtered data
genotypes <- read.table(file = "variants_maf5_miss9.012", header = FALSE, sep = "\t", 
                        na.strings = c("NA", "-1"))

# correct row names as sample id's
genotypes <- genotypes[, -1]
ind <- read.table(file = "variants_maf5_miss9.012.indv", header = FALSE)
rownames(genotypes) <- ind$V1

# read in parent offspring (F0 and Test F1 sample id's)
testsamp <- individuals$individual # 208 individuals

#create genotype matrix parent-offspring by filtering rownames of gmmat to include only those in testsamp$sample
genotypes <- genotypes[rownames(genotypes) %in% testsamp, , drop = FALSE]

# now that the dataframes have been called, there are formatting changes that need to be made.

# first, individuals$type needs to include the sex of the F0 parents
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
geno_conv <- as.data.frame(apply(genotypes, 2, function(x) {
  ifelse(x == 1, 2, ifelse(x == 2, 1, x))
}))

# the inspection done here lets you verify that all of your individuals are present
# in the genotype matrix
inspection <- inspect(ind = individuals, gen = geno_conv)
head(inspection)
inspection[which(inspection$sampled == 0), ] # got everyone!

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
geno_conv <- geno_conv[rownames(geno_conv) != "SAR_15_6433", ]
individuals <- individuals %>% 
    filter(individual != "SAR_15_6433")

# just need to remove it from three dataframes: geno_conv, individuals, and sar_clean. 207 obs -> inspection

# this is what generates for you all of the possible relationships at play
# this is to be subsetted downstream
combinations <- hothiphop(ind = individuals, gen = geno_conv)

# now we will return the top three most likely relationships for each offspring
# according to the hothiphop score
top3 <- topmatch(x = combinations, ranking = "hothiphop.parents")



## -- Exploration of the inferred relationships -- ##
## -- Rerun, without SAR_15_6433 -- ## 

# remove social parent entries (all NA anyway in this case)
top3 <- top3 %>% 
  filter(dam.type != "social parent")

# this neat chunk of code does two important things. 
# 1. if, for an individual, the Crossed_with entry is empty, turn it to an NA
# 2. take each Crossed_with entry, and if it consists of two individuals, split them
#    by the comma that separates the two sample id's. this prepares you for the
#    creation of the lookup table below
sar_clean <- sar_clean %>%
  mutate(Crossed_with_clean = ifelse(Crossed_with == "", NA, Crossed_with)) %>%
  mutate(Crossed_with_clean = strsplit(as.character(Crossed_with), ","))

# now we need to create a lookup table against which we will check all of the 
# mother-father relationships inferred by topmatch() and stored in top3
cross_lookup <- sar_clean %>%
  # keep only each sample, and who it's crossed with
  select(Sample_ID, Crossed_with_clean) %>%
  # this is fascinating here. this takes each of the instances where there are 3 
  # or more individuals included in a cross, and explodes them into their own row
  # of this lookup table. now instead of having a single row describe a cross between
  # 3 individuals, it expands to two (or three) rows of pairwise crosses
  unnest(Crossed_with_clean) %>%
  rename(partner = Crossed_with_clean) %>%
  filter(!is.na(partner)) %>%
  # this last line creates the pairwise cross, sorted so that the alphabetically 
  # smaller sample id appears first, and the two id's are separated by an underscore
  mutate(pair = paste(pmin(Sample_ID, partner), pmax(Sample_ID, partner), sep = "_"))

# some of the original entries in the sar_data dataframe included spaces.
# this code tidies that up in the lookup table using str_replace_all(spaces with no spaces)
# and does so in all columns.
# it ALSO creates a pair entry for each cross in the lookup table, where the cross is listed
# as being sorted where the smaller sample id is first, so it solves any ordering issues
library(stringr)
cross_lookup <- cross_lookup %>%
  mutate(
    Sample_ID = str_replace_all(str_trim(Sample_ID), " ", ""),
    partner = str_replace_all(str_trim(partner), " ", ""),
    pair = paste(pmin(Sample_ID, partner), pmax(Sample_ID, partner), sep = "_")
  )

# this is a ridiculously effective set of piped functions here:
# first, we create this new dataframe from top3
top3_checked <- top3 %>%
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
top3_checked <- top3_checked %>%
  select(year, brood, offspring, rank, dam, sire, pair, valid_cross, everything())

## -- Some summary information -- ##
length(unique(top3_checked$offspring)) # for 95 test F1's
table(top3_checked$valid_cross) # we're capturing 90 true relationships by using 
# the top 3 ranked pairs for each offspring

# among the top ranked parent pair for each offspring...
rank1 <- top3_checked %>% 
  filter(rank == 1)
table(rank1$valid_cross) # we captured the correct relationship 84/95 times 
# (88.42% accuracy)

# among the second highest ranked parent pair for each offspring...
rank2 <- top3_checked %>% 
  filter(rank == 2)
table(rank2$valid_cross) # we captured the correct relationship an additional 3 times...

# and among the third highest ranked parent pair for each offspring...
rank3 <- top3_checked %>% 
  filter(rank == 3) 
table(rank3$valid_cross) # we captured the correct relationship an additional 3 times..

# Alright, taking out SAR_15_6433 didn't affect our number of relationships.
# Here's the summary: (repeated from above)

# So for our 90/95 correct relationships...
# 84 came from the highest ranked pairs,
# 3 came from the second highest ranked pairs, and
# 3 came from the third highest ranked pairs.

# My question is... what happened to the other 5 offspring? Was it that:
# a) their true parent pairs are ranked lower than the third most likely?
# b) the hothiphop exclusion method didn't have the ability to recognize them at all?
# c) their parents were crossed but not documented? 
          # how were the crosses handled? was each cross incubated in its own container?
          # i would say that there's a reasonable chance that fertilization occurred after pooling


write.csv(top3_checked, file = "top3checked_f0testf1_043025.csv")

truefirst <- rank1 %>% 
    filter(valid_cross == TRUE)
summary(truefirst$hothiphop.parents)

falsefirst <- rank1 %>% 
  filter(valid_cross == FALSE)
summary(falsefirst$hothiphop.parents)

truesecond <- rank2 %>% 
  filter(valid_cross == TRUE)
summary(truesecond$hothiphop.parents)

falsesecond <- rank2 %>% 
  filter(valid_cross == FALSE)
summary(falsesecond$hothiphop.parents)

truethird <- rank3 %>% 
  filter(valid_cross == TRUE)
summary(truethird$hothiphop.parents)

falsethird <- rank3 %>% 
  filter(valid_cross == FALSE)
summary(falsethird$hothiphop.parents)

trueall <- top3_checked %>% 
    filter(valid_cross == TRUE)
summary(trueall$hothiphop.parents)

falseall <- top3_checked %>% 
  filter(valid_cross == FALSE)
summary(falseall$hothiphop.parents)


rank2true <- rank2 %>% 
    filter(valid_cross == TRUE)
rank3true <- rank3 %>% 
    filter(valid_cross == TRUE)
rank1false <- rank1 %>% 
    filter(valid_cross == FALSE)

tfpositions <- as.numeric(factor(top3_checked$valid_cross))
tfpositionsrank2true <- c(2, 2, 2)
tfpositionsrank3true <- c(2, 2, 2)
tfpositionsrank1false <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
boxplot(top3_checked$hothiphop.parents ~ top3_checked$valid_cross)
points(jitter(tfpositions, amount = 0.15), top3_checked$hothiphop.parents, pch = 19, col = tfpositions + 1)
points(jitter(tfpositionsrank2true, amount = 0.15), rank2true$hothiphop.parents, pch = 19, col = "blue")
points(jitter(tfpositionsrank3true, amount = 0.15), rank3true$hothiphop.parents, pch = 19, col = "blue")
points(jitter(tfpositionsrank1false, amount = 0.015), rank1false$hothiphop.parents, pch = 19, col = "blue")

# Okay, this tells me that these true relationships have high hothiphop scores. So another pair, even if
# it's still false, if it's got a lower score, it's going to be selected as the most likely relationship. 
# Let's see if filtering the data more strictly will affect the ability to detect the true relationship.
# For the top ranked relationships that are false, 

# I'm going to take the maf1 miss90 dataset, and I'm going to filter it for everything at once.
# (variants_maf1_miss9.recode.vcf is what generated the genotype matrix for this dataset.)

# Objective: Filter more strictly (with quality, and with quality + min/max depth) to see if
           # these additional filters remove sites that are causing problems with 















