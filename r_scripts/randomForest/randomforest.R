## randomforest.R by SPJ 062625
## PURPOSE: to run a randomforest model on SNPs aligned to YPE genome.
# goal is to see if we can predict sex from any of these SNPs (particularly Ch7)
# to have a chance to run hiphop on the F2 -> F1 assignment step
## USAGE: Rscript randomforest.R

# locally:
setwd("/Users/samjohnson/Desktop/")

# read in genotype matrix
mat <- read.table(file = "hard_variants_pflav_bial_mem_t2_noindels_q20_mindep3_maxdep75_maf1_miss80.012", header = FALSE, sep = "\t")

# correct row names as sample id's
mat<- mat[, -1]
gmmat <- as.matrix(mat)
ind <- read.table(file = "hard_variants_pflav_bial_mem_t2_noindels_q20_mindep3_maxdep75_maf1_miss80.012.indv", header = FALSE)

str(ind)
ind <- ind %>% 
  rename(sample = V1)
rownames(gmmat) <- ind$sample

pos <- read.table(file = "hard_variants_pflav_bial_mem_t2_noindels_q20_mindep3_maxdep75_maf1_miss80.012.pos", header = FALSE)
pos <- pos$V1
colnames(gmmat) <- pos

f0 <- read.csv(file = "spawned_f0.csv", header = FALSE)
colnames(f0) <- c("sample", "sex")

gmmat_f0 <- gmmat[rownames(gmmat) %in% f0$sample, ]

gmmat_f0 <- as.data.frame(gmmat_f0)
gmmat_f0 <- data.frame(sample = NA, gmmat_f0)
gmmat_f0$sample <- rownames(gmmat_f0)


gmmat_f0_sex <- cbind(gmmat_f0, sex = f0$sex)
gmmat_f0_sex <- gmmat_f0_sex %>% 
    select(sample, sex, everything())

install.packages("randomForest")
library(randomForest)

rf_model <- randomForest(as.factor(sex) ~ ., data = gmmat_f0_sex, importance = TRUE)
print(rf_model)
importance(rf_model)[,4]

without_sample <- gmmat_f0_sex %>% 
    select(sex, everything())

importance <- data.frame(importance(rf_model)[,4])
summary(importance)

