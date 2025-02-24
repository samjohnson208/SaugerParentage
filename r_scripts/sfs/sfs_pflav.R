getwd()
setwd("/Users/samjohnson/Documents/Sauger_102824/GeneticData/SFS/sfs_pflav/")
# read in data, remember to assign -1 as an NA string
var <- read.table(file = "variants_maf1_miss5.012", header = FALSE, sep = "\t", 
                  na.strings = c("NA", "-1"))

# remove the first column...
var <- var[, -1]

colnames(var)
dim(var) # notice that we've lost 9 sites... how'd that happen? tri/quadallelic sites?

sfs_pflav <- data.frame(site = 1:ncol(var), site_mean = NA, site_mean_d2 = NA, fold_site_mean_d2 = NA)
headtail(sfs_pflav)

# place site genotype means in sfs_svit$site_mean
for (i in 1:ncol(var)) {
  sfs_pflav$site_mean[i] <- mean(var[[i]], na.rm = TRUE)
}
# check
summary(sfs_pflav$site_mean)

# divide those by two to bound them [0,1], use these values for unfoldedSFS
for (i in 1:nrow(sfs_pflav)) {
  m <- sfs_pflav$site_mean[i]
  d2 <- m/2
  sfs_pflav$site_mean_d2[i] <- d2
}
# check
summary(sfs_pflav$site_mean_d2)

# if the mean value divided by two is greater than 0.5, subtract it from 1,
# and place that new value in the fold column, use for foldedSFS
for (i in 1:nrow(sfs_pflav)) {
  f <- sfs_pflav$site_mean_d2[i]
  if (f > 0.5){
    sfs_pflav$fold_site_mean_d2[i] <- 1-f
  } else {
    sfs_pflav$fold_site_mean_d2[i] <- f
  }
}

# check
summary(sfs_pflav$fold_site_mean_d2)
# this checks out! MAF should be 0.01 or higher, so the min is accurate!

# generate site frequency spectra 
# unfolded
hist(sfs_pflav$site_mean_d2, breaks = 100, main = "Unfolded SFS, YPE, MAF 0.01, Miss 0.5 (345 Sites)", xlab = "Allele Frequency")
# folded
hist(sfs_pflav$fold_site_mean_d2, breaks = 100, main = "Folded SFS, YPE, MAF 0.01, Miss 0.5 (345 Sites)", xlab = "Allele Frequency")




