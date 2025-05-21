### --- sauger_hiphop_allf1.R --- ###
### created by SPJ on 050925 ###
# purpose: explore hiphop output for ALL WILD F1'S AND ALL F0'S.
# compare that to the hiphop output for various PTPS scenarios for test indivs
# many more offspring included for this run, so we're pushing everything and running
# it on the cluster. all input files will go to...
# /project/ysctrout/hatchsauger/sam_sai_svit_mem/HipHop_Inp/HipHop_AllF1

# keep in mind, all of the setup to create these dataframes comes from the local script:
# /Users/samjohnson/Desktop/hiphop/sauger_hiphop_fifth.R

library(hiphop)
library(dplyr)
library(tidyverse)
library(purrr)

# read in the input files (these were scp'd to /project/ysctrout/hatchsauger/sam_sai_svit_mem/HipHop_Inp/HipHop_AllF1)
  individuals <- read.csv("wildF1_allF0_indivs.csv", header = TRUE, row.names = 1)
  hard_geno_conv <- read.csv("allF1F0genotypes.csv", header = TRUE, row.names = 1)

# run inspection
  f1_inspection <- inspect(ind = individuals, gen = hard_geno_conv)
  f1_inspection[which(f1_inspection$sampled == 0), ] # got everyone!
   
# now run combinations
  allf1_combinations <- hothiphop(ind = individuals, gen = hard_geno_conv)

# write out combinations to be processed locally
  write.csv(allf1_combinations, "allf1_combinations.csv")


