getwd()
setwd("/Users/samjohnson/Documents/Sauger_102824/Scripts_Plots/CountingReads/n_rawreads_per_indiv/")

reads <- read.table("nreads_per_ind.txt", header = FALSE, sep = ":")
colnames(reads) <- c("indiv", "nreads")
hist(reads$nreads, main = "number of raw reads per individual", xlab = "nreads per indiv")

reads <- data.frame(reads, thoureads = NA)
for(i in 1:nrow(reads)){
  n <- reads$nreads[i]
  reads$thoureads[i] <- n/1000
}
hist(reads$thoureads, main = "number of raw reads per individual", xlab = "nreads per indiv (thousands)")

reads <- data.frame(reads, milreads = NA)
for(i in 1:nrow(reads)){
  n <- reads$nreads[i]
  reads$milreads[i] <- n/1000000
}

hist(reads$milreads, main = "number of raw reads per individual (n = 1184)", xlab = "nreads/indiv (millions)")


