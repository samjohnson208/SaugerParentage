getwd()
setwd("/Users/samjohnson/Desktop/")

library(tidyverse)

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






##### ---- FROM WILL'S DATA. .FRQ INSTEAD OF .012 ---- #####
# read in data
wcr_data <- read.table(file = "rehead_variants_wcr_svit_mem_bial_noindels_q40_mindep8_maxdep75_miss95.frq", 
                       header = TRUE, sep = "\t")

# separate allele identities and frequencies
wcr_data <- wcr_data %>% 
    separate(X.ALLELE.FREQ1., into = c("allele1", "freq1"), sep = ":") %>% 
    separate(X.ALLELE.FREQ2., into = c("allele2", "freq2"), sep = ":")

# establish ref and alt frequencies
wcr_data$freq1 <- as.numeric(wcr_data$freq1)
wcr_data$freq2 <- as.numeric(wcr_data$freq2)

# hist(freq1) and hist(freq2) are opposite one another. we need an sfs for the minor allele
f2 <- wcr_data$freq2
wcr_data$maf <- pmin(f2, 1 - f2)
hist(wcr_data$maf)
hist(wcr_data$freq1)
hist(wcr_data$freq2)

# fit the beta to this sfs (josh's script)

#############################################
## empirical_beta_fit.R
#############################################

## JPJ 30 iv 24, modified by spj
## PURPOSE: to fit a beta distribution to a site frequency spectrum
## USAGE: Rscript empirical_beta_fit.R

install.packages("betafunctions")
library(betafunctions)

## read in data
#freqs <- read.delim("popN_125_1000_gen_20000_noHead.frq", header=FALSE)

## fold the allele frequencies
# folded <- matrix(NA, dim(freqs)[1], 1)
# for (i in 1:dim(freqs)[1]){
#   if (freqs[i,5] <= 0.5)	{ folded[i,1] <- freqs[i,5] }
#   else					{ folded[i,1] <- 1 - freqs[i,5] }
# }


# subset based on minor allele frequency cutoff
cutoff <- 0.0107
maf_cut <- matrix(NA, nrow(wcr_data), 1)
for (i in 1:nrow(wcr_data)) {
  if (wcr_data$maf[i] < cutoff)	{ maf_cut[i,1] <- 0 }
  else					{ maf_cut[i,1] <- 1 }
}
wcr_data$maf_include <- maf_cut

# include only loci that have maf over 0.01
wcr_data_maf_cut_folded <- wcr_data %>% 
    filter(maf_cut == 1)

# NOT DOING THIS NOW... RIGHT? THEY WANTED ME TO AVOID MESSING W MAF. YOU'RE THEORETICALLLY REMOVING
# MOST OF THE GENOTYPING ERROR BY REQUIRING A HIGH DEPTH... RIGHT?

## fit the 4 parameter beta (alpha, beta, lower, upper)
# beta_fit_maf_cut_folded <- Beta.4p.fit(maf_cut_folded[,1])
beta_fit_maf_cut_folded <- Beta.4p.fit(wcr_data_maf_cut_folded$maf)

# this is the key. this is what we need. i'm almost sure of it.
beta_fit_alt_freq <- Beta.4p.fit(wcr_data_maf_cut_folded$freq2)

## plot SFS histogram with fitted beta
# pdf("beta_fit_popN_125_1000_gen_20000.pdf", height=5, width=5)
# par(mar=c(5,5,1,1))
hist(wcr_data_maf_cut_folded$maf, breaks=100, xlab="Folded allele frequency", ylab="Number of loci", main="", cex.axis=1.25, cex.lab=1.5); box(lwd=1.5)
par(new=TRUE)
dbeta_seq_x <- seq(beta_fit_maf_cut_folded$l, beta_fit_maf_cut_folded$u, length=100)
dbeta_seq_y <- dBeta.4P(seq(beta_fit_maf_cut_folded$l, beta_fit_maf_cut_folded$u, length=100), beta_fit_maf_cut_folded$l, beta_fit_maf_cut_folded$u, beta_fit_maf_cut_folded$alpha, beta_fit_maf_cut_folded$beta)
plot(dbeta_seq_x[2:99], dbeta_seq_y[2:99], ylim=c(0,dbeta_seq_y[2]), type="l", col="red", lwd=2, xlab="", ylab="", axes=FALSE)
mtext(bquote("maf"==.(cutoff)), cex=1.25, adj=0.98, line=-1.25)
mtext(bquote(alpha==.(round(beta_fit_maf_folded$alpha,3))), cex=1.25, adj=0.98, line=-2.4)
mtext(bquote(beta==.(round(beta_fit_maf_folded$beta,3))), cex=1.25, adj=0.98, line=-3.75)
dev.off()

# alright, so what happened here was that the beta fit without that maf cutoff 
# looked just about horrible, so i added it. i still don't think it looks great,
# see how we're underestimating the sfs at those lower frequencies, but we're
# overestimating at those intermediate frequencies... gonna go ahead with the params
# estimated from the sfs with the 0.0107 cutoff for now. we'll discuss with them.








