## -- 02/18/25 -- ##

# This script was used on the cluster to generate histograms of the percentage 
# of missing data per site for the genetic datasets that were aligned to walleye 
# and yellow perch reference genomes, respectively.

# Notice the 1. higher number of reads for walleye as well as the 2. higher
# proportion of reads with no missing data in the walleye dataset. This comparison
# was made simply by visually analyzing the height of the histogram bars as compared
# to the y-axis ticks.

setwd("/project/ysctrout/hatchsauger/sam_sai_pflav/")
pflav <- read.table(file = "/project/ysctrout/hatchsauger/sam_sai_pflav/per_md_per_site_pflav.missing.lmiss", sep = "\t", header = TRUE)
head(pflav)
hist(pflav$F_MISS, breaks = 20, main = "% Missing Data per Site (Yellow Perch)", xlab = "% Missing Data per Site", ylim = c(0,200000))

setwd("/project/ysctrout/hatchsauger/sam_sai_svit/")
svit <- read.table(file = "/project/ysctrout/hatchsauger/sam_sai_svit/per_md_per_site_svit.missing.lmiss", sep = "\t", header = TRUE)
head(svit)
hist(svit$F_MISS, breaks = 20, main = "% Missing Data per Site (Walleye)", xlab = "% Missing Data per Site")