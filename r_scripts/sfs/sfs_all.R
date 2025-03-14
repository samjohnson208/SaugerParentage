getwd()
setwd("/Users/samjohnson/Desktop/")
# read in data, remember to assign -1 as an NA string
svitaln <- read.table(file = "variants_maf1_miss9_aln_svit.012_conv", header = FALSE, sep = "\t", 
                  na.strings = c("NA", "-9"))
pflavaln <- read.table(file = "variants_maf1_miss9_aln_pflav.012_conv", header = FALSE, sep = "\t", 
                  na.strings = c("NA", "-9"))
svitmem <- read.table(file = "variants_maf1_miss9_mem_svit.012_conv", header = FALSE, sep = "\t", 
                  na.strings = c("NA", "-9"))
pflavmem <- read.table(file = "variants_maf1_miss9_mem_pflav.012_conv", header = FALSE, sep = "\t", 
                  na.strings = c("NA", "-9"))

# remove the first column...
svitaln<- svitaln[, -1]
pflavaln<- pflavaln[, -1]
svitmem <- svitmem[, -1]
pflavmem <- pflavmem[, -1]

sfs_svit_aln <- data.frame(site = 1:ncol(svitaln), site_mean = NA, site_mean_d2 = NA, fold_site_mean_d2 = NA)
sfs_pflav_aln <- data.frame(site = 1:ncol(pflavaln), site_mean = NA, site_mean_d2 = NA, fold_site_mean_d2 = NA)
sfs_svit_mem <- data.frame(site = 1:ncol(svitmem), site_mean = NA, site_mean_d2 = NA, fold_site_mean_d2 = NA)
sfs_pflav_mem <- data.frame(site = 1:ncol(pflavmem), site_mean = NA, site_mean_d2 = NA, fold_site_mean_d2 = NA)

# place site genotype means in df$site_mean
for (i in 1:ncol(svitaln)) {
  sfs_svit_aln$site_mean[i] <- mean(svitaln[[i]], na.rm = TRUE)
}
for (i in 1:ncol(pflavaln)) {
  sfs_pflav_aln$site_mean[i] <- mean(pflavaln[[i]], na.rm = TRUE)
}
for (i in 1:ncol(svitmem)) {
  sfs_svit_mem$site_mean[i] <- mean(svitmem[[i]], na.rm = TRUE)
}
for (i in 1:ncol(pflavmem)) {
  sfs_pflav_mem$site_mean[i] <- mean(pflavmem[[i]], na.rm = TRUE)
}

# check
summary(sfs_svit_aln$site_mean)
summary(sfs_pflav_aln$site_mean)
summary(sfs_svit_mem$site_mean)
summary(sfs_pflav_mem$site_mean)


# divide those by two to bound them [0,1], use these values for unfoldedSFS
for (i in 1:nrow(sfs_svit_aln)) {
  m <- sfs_svit_aln$site_mean[i]
  d2 <- m/2
  sfs_svit_aln$site_mean_d2[i] <- d2
}
for (i in 1:nrow(sfs_pflav_aln)) {
  m <- sfs_pflav_aln$site_mean[i]
  d2 <- m/2
  sfs_pflav_aln$site_mean_d2[i] <- d2
}
for (i in 1:nrow(sfs_svit_mem)) {
  m <- sfs_svit_mem$site_mean[i]
  d2 <- m/2
  sfs_svit_mem$site_mean_d2[i] <- d2
}
for (i in 1:nrow(sfs_pflav_mem)) {
  m <- sfs_pflav_mem$site_mean[i]
  d2 <- m/2
  sfs_pflav_mem$site_mean_d2[i] <- d2
}

# check
summary(sfs_svit_aln$site_mean_d2)
summary(sfs_pflav_aln$site_mean_d2)
summary(sfs_svit_mem$site_mean_d2)
summary(sfs_pflav_mem$site_mean_d2)

# if the mean value divided by two is greater than 0.5, subtract it from 1,
# and place that new value in the fold column, use for foldedSFS
for (i in 1:nrow(sfs_svit_aln)) {
  f <- sfs_svit_aln$site_mean_d2[i]
  if (f > 0.5){
    sfs_svit_aln$fold_site_mean_d2[i] <- 1-f
  } else {
    sfs_svit_aln$fold_site_mean_d2[i] <- f
  }
}

for (i in 1:nrow(sfs_pflav_aln)) {
  f <- sfs_pflav_aln$site_mean_d2[i]
  if (f > 0.5){
    sfs_pflav_aln$fold_site_mean_d2[i] <- 1-f
  } else {
    sfs_pflav_aln$fold_site_mean_d2[i] <- f
  }
}

for (i in 1:nrow(sfs_svit_mem)) {
  f <- sfs_svit_mem$site_mean_d2[i]
  if (f > 0.5){
    sfs_svit_mem$fold_site_mean_d2[i] <- 1-f
  } else {
    sfs_svit_mem$fold_site_mean_d2[i] <- f
  }
}

for (i in 1:nrow(sfs_pflav_mem)) {
  f <- sfs_pflav_mem$site_mean_d2[i]
  if (f > 0.5){
    sfs_pflav_mem$fold_site_mean_d2[i] <- 1-f
  } else {
    sfs_pflav_mem$fold_site_mean_d2[i] <- f
  }
}

# check
summary(sfs_svit_aln$fold_site_mean_d2)
summary(sfs_pflav_aln$fold_site_mean_d2)
summary(sfs_svit_mem$fold_site_mean_d2)
summary(sfs_pflav_mem$fold_site_mean_d2)

# this checks out! MAF should be 0.01 or higher, so the min is accurate!

# generate site frequency spectra, unfolded
par(mfrow=c(2,2))
hist(sfs_svit_aln$site_mean_d2, breaks = 100, main = "Unfolded SFS, WAE ALN (6468 Sites)", xlab = "Ref Allele Frequency")
hist(sfs_pflav_aln$site_mean_d2, breaks = 100, main = "Unfolded SFS, YPE ALN (166 Sites)", xlab = "Ref Allele Frequency")
hist(sfs_svit_mem$site_mean_d2, breaks = 100, main = "Unfolded SFS, WAE MEM (7369 Sites)", xlab = "Ref Allele Frequency")
hist(sfs_pflav_mem$site_mean_d2, breaks = 100, main = "Unfolded SFS, YPE MEM (204 Sites)", xlab = "Ref Allele Frequency")

# generate site frequency spectra, folded
par(mfrow=c(2,2))
hist(sfs_svit_aln$fold_site_mean_d2, breaks = 100, main = "Folded SFS, WAE ALN (6468 Sites)", xlab = "Ref Allele Frequency")
hist(sfs_pflav_aln$fold_site_mean_d2, breaks = 100, main = "Folded SFS, YPE ALN (166 Sites)", xlab = "Ref Allele Frequency")
hist(sfs_svit_mem$fold_site_mean_d2, breaks = 100, main = "Folded SFS, WAE MEM (7369 Sites)", xlab = "Ref Allele Frequency")
hist(sfs_pflav_mem$fold_site_mean_d2, breaks = 100, main = "Folded SFS, YPE MEM (204 Sites)", xlab = "Ref Allele Frequency")

# # unfolded
# hist(sfs_svit$site_mean_d2, breaks = 100, main = "Unfolded SFS, WAE, MAF 0.01, Miss 0.9 (6468 Sites)", xlab = "Allele Frequency")
# # folded
# hist(sfs_svit$fold_site_mean_d2, breaks = 100, main = "Folded SFS, WAE, MAF 0.01, Miss 0.9 (6468 Sites)", xlab = "Allele Frequency")



