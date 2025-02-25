install.packages("sequoia")
library(sequoia)
library(dplyr)
?sequoia()

getwd()
setwd("/Users/samjohnson/Documents/Sauger_102824/GeneticData/Sequoia_Inp")

mat <- read.table(file = "variants_maf1_miss9.012_conv", header = FALSE, sep = "\t", 
                  na.strings = c("NA", "-1"))
mat<- mat[, -1]
gmmat <- as.matrix(mat)

ind <- read.table(file = "variants_maf1_miss9.012.indv", header = FALSE)
str(ind)
ind <- ind %>% 
    rename(sample = V1)

rownames(gmmat) <- ind$sample

check <- CheckGeno(gmmat, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))    
# ✖ There are 2 individuals scored for <5% of SNPs, these WILL BE IGNORED
# ✖ In addition, there are 6 individuals scored for <20% of SNPs,it is advised to treat their assignments with caution
# ℹ After exclusion, There are  1182  individuals and  6468  SNPs.

stats <- SnpStats(gmmat, Pedigree = NULL, Duplicates = NULL, Plot = TRUE, quiet = FALSE)
# ℹ Not conditioning on any pedigree



