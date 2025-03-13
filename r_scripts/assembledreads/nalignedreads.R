getwd()
setwd("/Users/samjohnson/Documents/Sauger_102824/GeneticData/AssembledReads/n_aln_per_indiv/")

## Walleye Reference Genome ##
alnwall <- read.table("assembled_per_ind_svit.txt", header = TRUE, sep = " ")
head(alnwall)
dim(alnwall)

#in numbers of reads
hist(alnwall$assembled, main = "number of aligned reads per individual (Walleye Reference)", xlab = "nreads/indiv")

#create a column for thousands of aligned reads
alnwall <- data.frame(alnwall, thou_aln_reads = NA)
for(i in 1:nrow(alnwall)){
  n <- alnwall$assembled[i]
  alnwall$thou_aln_reads[i] <- n/1000
}

# in thousands of reads
hist(alnwall$thou_aln_reads, main = "number of aligned reads per individual (Walleye Reference)", xlab = "nreads per indiv (thousands)")


## Yellow Perch Reference Genome ##
alnype <- read.table("assembled_per_ind_pflav.txt", header = TRUE, sep = " ")
head(alnype)
dim(alnype)

#in numbers of reads
hist(alnype$assembled, main = "number of aligned reads per individual (Yellow Perch Reference)", xlab = "nreads/indiv")

#create a column for thousands of aligned reads
alnype <- data.frame(alnype, thou_aln_reads = NA)
for(i in 1:nrow(alnype)){
  n <- alnype$assembled[i]
  alnype$thou_aln_reads[i] <- n/1000
}

hist(alnype$thou_aln_reads, main = "number of aligned reads per individual (Yellow Perch Reference)", xlab = "nreads per indiv (thousands)")

par(mfrow=c(1,1))
hist(alnwall$thou_aln_reads, main = "number of aligned reads per individual (Walleye and Yellow Perch Reference Genomes)", 
     xlab = "nreads per indiv (thousands)", col = "darkgreen", breaks = 15, ylim = c(0,350))
hist(alnype$thou_aln_reads, col = "orange", add = TRUE, breaks = 15)
legend("topright", legend=c("Walleye", "Yellow Perch"), fill=c("darkgreen", "orange"))


#now for % aligned
alnwall <- data.frame(alnwall, per_aln = NA)
alnype <- data.frame(alnype, per_aln = NA)

for(i in 1:nrow(alnwall)){
  nassem <- alnwall$assembled[i]
  nraw <- alnwall$raw[i]
  alnwall$per_aln[i] <- nassem/nraw*100
}

for(i in 1:nrow(alnype)){
  nassem <- alnype$assembled[i]
  nraw <- alnype$raw[i]
  alnype$per_aln[i] <- nassem/nraw*100
}
par(mfrow = c(1,1))
hist(alnwall$per_aln, main = "percent raw reads aligned to the reference by individual (n = 1184)", 
     xlab = "%rawreadsaligned/individual", col = "darkgreen", breaks = 15, ylim = c(0,800), xlim = c(0,60))
hist(alnype$per_aln, col = "orange", add = TRUE, breaks = 15)
legend("topright", legend=c("Walleye", "Yellow Perch"), fill=c("darkgreen", "orange"))

summary(alnype$per_aln)

# now to add european perch (pfluv)

getwd()
setwd("/Users/samjohnson/Documents/Sauger_102824/GeneticData/AssembledReads/n_aln_per_indiv/")

euro <- read.table("assembled_per_ind_pfluv.txt", header = TRUE, sep = " ")

hist(euro$assembled, main = "number of aligned reads per individual (Euro Perch Reference)", xlab = "nreads/indiv")

euro <- data.frame(euro, thou_aln_reads = NA)
for(i in 1:nrow(euro)){
  n <- euro$assembled[i]
  euro$thou_aln_reads[i] <- n/1000
}

hist(euro$thou_aln_reads, main = "number of aligned reads per individual (Walleye Reference)", xlab = "nreads per indiv (thousands)")

hist(alnwall$thou_aln_reads, main = "number of aligned reads per individual by reference genome", 
     xlab = "nreads per indiv (thousands)", col = "darkgreen", breaks = 15, ylim = c(0,350))
hist(alnype$thou_aln_reads, col = "orange", add = TRUE, breaks = 15)
hist(euro$thou_aln_reads, col = "skyblue", add = TRUE, breaks = 15)
legend("topright", legend=c("Walleye", "Yellow Perch", "European Perch"), fill=c("darkgreen", "orange", "skyblue"))

euro <- data.frame(euro, per_aln = NA)
for(i in 1:nrow(euro)){
  nassem <- euro$assembled[i]
  nraw <- euro$raw[i]
  euro$per_aln[i] <- nassem/nraw*100
}

hist(alnwall$per_aln, main = "percent raw reads aligned per individual by reference genome (n = 1184)", 
     xlab = "%rawreadsaligned/individual", col = "darkgreen", breaks = 15, ylim = c(0,800), xlim = c(0,60))
hist(alnype$per_aln, col = "orange", add = TRUE, breaks = 15)
hist(euro$per_aln, col = "skyblue", add = TRUE, breaks = 15)
legend("topleft", legend=c("Walleye", "Yellow Perch", "European Perch"), fill=c("darkgreen", "orange", "skyblue"))


# checking out WR pflav

getwd()
setwd("/Users/samjohnson/Documents/Sauger_102824/GeneticData/AssembledReads/n_aln_per_indiv/")

## WR YPE Reference Genome ##
wrype <- read.table("assembled_per_ind_pflav_WR.txt", header = TRUE, sep = " ")
head(wrype)

#create a column for thousands of aligned reads
wrype <- data.frame(wrype, thou_aln_reads = NA)

for(i in 1:nrow(wrype)){
  n <- wrype$assembled[i]
  wrype$thou_aln_reads[i] <- n/1000
}

# now for % aligned
wrype <- data.frame(wrype, per_aln = NA)

for(i in 1:nrow(wrype)){
  nassem <- wrype$assembled[i]
  nraw <- wrype$raw[i]
  wrype$per_aln[i] <- nassem/nraw*100
}


hist(wrype$per_aln, main = "percent raw reads aligned to the reference by individual (n = 1184)", 
     xlab = "%rawreadsaligned/individual", col = "salmon", breaks = 15)
hist(alnype$per_aln, col = "skyblue", add = TRUE, breaks = 15)
legend("topleft", legend=c("WR - YPE", "SJ - YPE"), fill=c("salmon", "skyblue"))


# checking out BWA MEM output
getwd()
setwd("/Users/samjohnson/Documents/Sauger_102824/GeneticData/AssembledReads/n_aln_per_indiv/")

wae <- read.table("assembled_per_ind_svit.txt", header = TRUE, sep = " ")
ype <- read.table("assembled_per_ind_pflav.txt", header = TRUE, sep = " ")

wae_mem <- read.table("assembled_per_ind_svit_mem.txt", header = TRUE, sep = " ")
ype_mem <- read.table("assembled_per_ind_pflav_mem.txt", header = TRUE, sep = " ")

# add column for percent aligned
wae <- data.frame(wae, per_aln = NA)
ype <- data.frame(ype, per_aln = NA)
wae_mem <- data.frame(wae_mem, per_aln = NA)
ype_mem <- data.frame(ype_mem, per_aln = NA)

for(i in 1:nrow(wae)){
  nassem <- wae$assembled[i]
  nraw <- wae$raw[i]
  wae$per_aln[i] <- nassem/nraw*100
}

for(i in 1:nrow(ype)){
  nassem <- ype$assembled[i]
  nraw <- ype$raw[i]
  ype$per_aln[i] <- nassem/nraw*100
}

for(i in 1:nrow(wae_mem)){
  nassem <- wae_mem$assembled[i]
  nraw <- wae_mem$raw[i]
  wae_mem$per_aln[i] <- nassem/nraw*100
}

for(i in 1:nrow(ype_mem)){
  nassem <- ype_mem$assembled[i]
  nraw <- ype_mem$raw[i]
  ype_mem$per_aln[i] <- nassem/nraw*100
}


hist(wae_mem$per_aln, col = "springgreen2", main = "percent raw reads aligned to the reference by individual (n = 1184)", 
     xlab = "% rawreads aligned", breaks = 25)
hist(wae$per_aln, col = "darkgreen", add = TRUE, breaks = 25)
hist(ype$per_aln, col = "goldenrod2", add = TRUE, breaks = 25)
hist(ype_mem$per_aln, col = "yellow2", add = TRUE, breaks = 25)
legend("topleft", legend=c("WAE - aln/samse", "WAE - mem", 
                           "YPE - aln/samse", "YPE - mem"), 
       fill=c("darkgreen", "springgreen2", "goldenrod2", "yellow2"))






