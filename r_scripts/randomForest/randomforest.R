## randomforest.R by SPJ 062625
## PURPOSE: to run a randomforest model on SNPs aligned to YPE genome.
# goal is to see if we can predict sex from any of these SNPs (particularly Ch7)
# to have a chance to run hiphop on the F2 -> F1 assignment step
## USAGE: Rscript randomforest.R

# packages
install.packages("randomForest")
library(randomForest)
library(tidyverse)

# locally:
setwd("/Users/samjohnson/Desktop/")

# read in genotype matrix
mat <- read.table(file = "hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss95.012", 
                  header = FALSE, sep = "\t")

# correct row names as sample id's
mat<- mat[, -1]
gmmat <- as.matrix(mat)

# read in individuals for rownames(gmmat)
ind <- read.table(file = "hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss95.012.indv", 
                  header = FALSE)
str(ind)
ind <- ind %>% 
  rename(sample = V1)
rownames(gmmat) <- ind$sample # add inds to gmmat

# read in positions for colnames(gmmat)
pos <- read.table(file = "hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss95.012.pos", 
                  header = FALSE)
pos$V3 <- paste(pos$V1, pos$V2, sep = "_") # combine scaffold and position to one entry
colnames(gmmat) <- pos$V3

# read in list of F0 individuals and their sexes
f0 <- read.csv(file = "spawned_f0.csv", header = FALSE)
colnames(f0) <- c("sample", "sex")

# filter gmmat to include only F0's
gmmat_f0 <- gmmat[rownames(gmmat) %in% f0$sample, ]
gmmat_f0 <- as.data.frame(gmmat_f0) # convert to df

# need to add this sample column to cbind, going to have to delete after
gmmat_f0 <- data.frame(sample = NA, gmmat_f0)
gmmat_f0$sample <- rownames(gmmat_f0)

# can we cbind()? are they in the same order?
table(f0$sample == rownames(gmmat_f0))

# cbind() and check
gmmat_f0_sex <- cbind(gmmat_f0, sex = f0$sex)
table(f0$sex == gmmat_f0_sex$sex) # good!

# bring sex to the front of the df and remove sample so it is not a part of the analysis
gmmat_f0_sex <- gmmat_f0_sex %>% 
    select(sample, sex, everything())
gmmat_f0_sex <- gmmat_f0_sex[, -1]

# construct randomForest model
rf_model <- randomForest(as.factor(sex) ~ ., data = gmmat_f0_sex, importance = TRUE)
print(rf_model)
importance(rf_model)[,4]

importance <- data.frame(importance(rf_model)[,4])
summary(importance)
# unsuccessful.
