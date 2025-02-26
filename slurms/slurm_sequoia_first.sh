#!/bin/bash

## created by SPJ 022525
## PURPOSE: run sequoia on R on the cluster (on a converted genotype matrix)
## USAGE: sbatch slurm_sequoia_first.sh

#SBATCH --job-name=sequoia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=10:00:00
#SBATCH --mem=10G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_sequoia_first
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_sequoia_first

# cd to directory with your converted genotype matrix
cd /project/ysctrout/hatchsauger/sam_sai_svit/Sequoia_Inp

# load modules
module load arcc/1.0 gcc/14.2.0 r/4.4.0
R --save << EOF

# unsure if this will work first try. worth a shot.
install.packages("sequoia")
65
library(sequoia)

install.packages("dplyr")
65
qlibrary(dplyr)

sink("sequoia_cluster_out.R")

# read in genotype matrix
mat <- read.table(file = "variants_maf1_miss9.012_conv", header = FALSE, sep = "\t", 
                  na.strings = c("NA", "-1"))

# correct row names as sample id's
mat<- mat[, -1]
gmmat <- as.matrix(mat)
ind <- read.table(file = "variants_maf1_miss9.012.indv", header = FALSE)
str(ind)
ind <- ind %>% 
    rename(sample = V1)
rownames(gmmat) <- ind$sample

# check genotype matrix for samples/loci to be excluded
check <- CheckGeno(gmmat, quiet = FALSE, Plot = TRUE, Return = "GenoM", Strict = TRUE, DumPrefix = c("F0", "M0"))    
# There are 2 individuals scored for <5% of SNPs, these WILL BE IGNORED
# In addition, there are 6 individuals scored for <20% of SNPs,it is advised to treat their assignments with caution
# After exclusion, There are  1182  individuals and  6468  SNPs.

# generate snp stats and plots for the genotype matrix
stats <- SnpStats(gmmat, Pedigree = NULL, Duplicates = NULL, Plot = TRUE, quiet = FALSE)
# Not conditioning on any pedigree

# run sequoia on the genotype matrix
outfull <- sequoia(GenoM = gmmat, Module = 'ped', MaxSibIter = 42, StrictGenoCheck = TRUE, CalcLLR = TRUE)

# run GetMaybeRel() on the sequoia output
gmr <- GetMaybeRel(outfull, GenoM = gmmat)


# save output to current wd 
save(outfull, file = "Sequoia_OutFull_022525.RData")
save(gmr, file = "Sequoia_GetMayRel_022525.RData")


sink()


EOF

#exit R
q()





