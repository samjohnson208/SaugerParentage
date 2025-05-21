### --- sauger_hiphop_fifth.R --- ###
### created by SPJ on 050925 ###
# purpose: explore hiphop output for ALL WILD F1'S AND ALL F0'S.
# compare that to the hiphop output for various PTPS scenarios on test indivs

getwd()
setwd("/Users/samjohnson/Desktop/hiphop")

#install.packages("hiphop")
library(hiphop)
library(dplyr)
library(tidyverse)
library(purrr)

# read in test results
top1_100_75_50_25 <- read.csv("top1_100_75_50_25.csv", header = TRUE)
# observations for 85 of the test offspring at varying PTPS

# read in individuals .csv
individuals <- read.csv("wildF1_allF0_indivs.csv", header = TRUE, row.names = 1)

# read in genotypes
hard_genotypes <- read.table(file = "hard_variants_bial_noindels_q20_mindep8_maxdep75_maf5_miss95.012", 
                             header = FALSE, sep = "\t", na.strings = c("NA", "-1"))
hard_genotypes <- hard_genotypes[, -1]
hard_ind <- read.table(file = "hard_variants_bial_noindels_q20_mindep8_maxdep75_maf5_miss95.012.indv", header = FALSE)
rownames(hard_genotypes) <- hard_ind$V1

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

hard_geno_conv <- as.data.frame(apply(hard_genotypes, 2, function(x) {
  ifelse(x == 1, 2, ifelse(x == 2, 1, x))
}))

# remove that one parent that wasn't spawned
individuals <- individuals %>% 
  filter(individual != "SAR_15_6433")
hard_geno_conv <- hard_geno_conv[rownames(hard_geno_conv) != "SAR_15_6433", ]

# create a vector of all individuals to be included in this analysis
indivs_wildf1_allf0 <- individuals$individual # 596 individuals

# create genotype matrix for only these parents/offspring by filtering rownames of 
# gmmat to include only those in indivs_wildf1_allf0
hard_geno_conv <- hard_geno_conv[rownames(hard_geno_conv) %in% indivs_wildf1_allf0, , drop = FALSE]

length(indivs_wildf1_allf0)
dim(hard_geno_conv)
non_matching_ids <- !indivs_wildf1_allf0 %in% rownames(hard_geno_conv)
indivs_wildf1_allf0[non_matching_ids]
# > indivs_wildf1_allf0[non_matching_ids]
# [1] "SAR_21_5582" "SAR_21_5641" "SAR_21_5647" "SAR_21_5652" "SAR_21_5771" "SAR_21_5819"
# this is probably because i pulled from sar_data rather than extractions.
# no 5582 fin clip, same for 5641, 5647, 5652, 5771, no vial for 5819
# i mean, this is exactly what the inspection step will do...
# let's filter those out of indivs_wildf1_allf0

individuals <- individuals %>% 
    filter(individual %in% indivs_wildf1_allf0)
table(individuals$type)
# nice. 590 obs. 328 f1 offspring (juvenile AND spawning) 262 parents (2015 and 2016) 

# cleaning is now complete. we're ready to start hiphop functions!
f1_inspection <- inspect(ind = individuals, gen = hard_geno_conv)
f1_inspection[which(f1_inspection$sampled == 0), ] # got everyone!

write.csv(individuals, "wildF1_allF0_indivs.csv")
write.csv(hard_geno_conv, "allF1F0genotypes.csv")

# To run on the cluster!
# individuals <- read.csv("wildF1_allF0_indivs.csv", header = TRUE)
# hard_geno_conv <- read.csv("allF1F0genotypes.csv", header = TRUE)

# run hothiphop
allf1_combinations <- hothiphop(ind = individuals, gen = hard_geno_conv)
# five MILLION observations. wow.
write.csv(allf1_combinations, "allf1_combinations.csv")

# keep the top 3 parent relationships for each offspring, remove social parent rows
top_f1 <- topmatch(x = allf1_combinations, ranking = "hothiphop.parents")
top_f1 <- top_f1 %>% 
  filter(dam.type != "social parent")

# filter for only top ranked relationships
top_f1 <- top_f1 %>%
  filter(rank == 1)

# now generates the inferred pair column for this ptps
top_f1 <- top_f1 %>% 
  mutate(pair = paste(pmin(dam, sire), pmax(dam, sire), sep = "_"))

write.csv(top_f1, "top_f1.csv")


