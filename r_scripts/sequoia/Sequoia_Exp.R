install.packages("sequoia")
library(sequoia)
library(dplyr)


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

# now let's take a look at the MAFs and attempt to construct a dataframe with well
# filtered data that will maximize computational efficiency
stats <- as.data.frame(stats)
stats <- data.frame(stats$AF, MAF = NA, stats$Mis, stats$HWE.p)
for (i in 1:nrow(stats)) {
  f <- stats$stats.AF[i]
  if (f > 0.5){
    stats$MAF[i] <- 1-f
  } else {
    stats$MAF[i] <- f
  }
}

# filter for snps that have higher than 20% MAF
statsmaf20 <- stats %>%
    filter(MAF > 0.2)
# 2228 loci remaining.

# create a reduced matrix with only 100 rows and columns
gmmat_red <- gmmat[1:100,1:100]

# run sequoia on the reduced matrix
sequoia_red <- sequoia(GenoM = gmmat_red, Module = 'ped', MaxSibIter = 42, StrictGenoCheck = TRUE, CalcLLR = TRUE)
save(SeqList, LHdata, Geno, file = "Sequoia_RedOut_022525.RData")

# Notes from first run output (reduced matrix)
  # Obviously not many relationships because these are Will's/Liz's fish and some of mine from 1 gen.

# See specs stored in the sequoia_red object as well as ErrM (genotyping error matrix)
# Still not understanding args.AP
# Note excluded monomorphic snps
# Need to understand format for AgePriors
# May prove extremely beneficial to add sex and age data. (see sequoia_red$LifeHist)

# Pedigree dataframe: some dummy parents assigned. Also stored in $DummyIDs, maybe some sibs there?
# May need to revisit the paper to refresh on the TotLikPar/Sib objects.

# Summary: Lots going on here, some of which I get and some of which I don't. 
# May help to see it run on the whole dataset. See how the total likelihood objects
# change. However, the main thing I want to see is whether or not more relationships are
# going to be formed by adding these other generations in there, even without the LH (sex, b.y.) data.

# Also need to understand what's going on with the plots and how they're generated.


GetMaybeRel(sequoia_red, GenoM = gmmat_red)

# Important to note: we've got a hundred loci here, so the OH across only a hundred loci
# is going to create a lot more potential relationships I think than are feasible. Yet
# another reason why I want to see this on the whole dataset.







