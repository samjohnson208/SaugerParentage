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



