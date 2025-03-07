getwd()
setwd("/Users/samjohnson/Documents/Sauger_102824/GeneticData/RadMapReport/RadMap_svit/") 

# The slurm_radmapreport_fix.sh outputs a .txt file that has 7 columns of interest:
# Sample, Library, Assembled Reads, lociN, meanDepth, lociNmin10reads, meanDepthMin10reads
# However, only 5 of these columns are filled with info (not Library and not Assembled Reads)

# Fortunately, I've calculated the number of assembled reads per individual already, using count_assembled.sh
# from Will. Before using this script, I've taken the number of assembled reads per individual from here:

# /Users/samjohnson/Documents/Sauger_102824/GeneticData/AssembledReads/n_aln_per_indiv

# made sure that sample id columns matched, and pasted those raw reads and assembled reads into empty columns of
# the radmap report. (excel) That generated mappingreport_svit_reads.txt. However, we still need to add the library
# for each sample to test for library effects in the analysis of nreads and read depth...

# There exists a .txt file that was extracted from the GTL form called sample_sampleLib.txt.
# That file has been placed in both the RadMap_svit and RadMap_pflav directories. (See path below)
# This contains the sample id and the plate/library that the sample was sequenced in.

# /Users/samjohnson/Documents/Sauger_102824/GeneticData/RadMapReport

##### ---------- Walleye Alignment ---------- #####

library(dplyr)

# load in both tables
svit <- read.table(file = "mappingreport_svit_reads.txt", header = TRUE, sep = "\t", na.strings = "")
samp <- read.table(file = "sample_sampleLib.txt", header = TRUE, sep = "\t")

# merge by sample id, conserve all rows
radmap_svit <- merge(samp, svit, by = "sample", all = TRUE)

# the sampleLib column from samp contains both the plate and library info.
radmap_svit <- radmap_svit %>% 
    rename(plate = sampleLib)

# create a new column to JUST store the library, so that we can color by lib in downstream plots
radmap_svit <- data.frame(radmap_svit, lib = NA)

# rearrange columns so that plate and library info are at the left of the df
radmap_svit <- radmap_svit[, c("sample", "plate", "lib", "rawreads", "assembledreads", "lociN", "meanDepth", "lociNmin10reads", "meanDepthMin10reads")]

# for every row in the df, if the plate is from lib 1, put 1 in the lib column, else 2 for lib 2.
for(i in 1:nrow(radmap_svit)) {
  if(grepl("Lib1", radmap_svit$plate[i])) {
    radmap_svit$lib[i] <- 1
  } else {
    radmap_svit$lib[i] <- 2
  }
}

# Now the full dataframe has been established as we wish it to be.
write.table(radmap_svit, file = "mappingreport_svit_reads_lib.txt")

plot(radmap_svit$lociNmin10reads ~ radmap_svit$assembledreads, col = radmap_svit$lib,
     main = "Loci with 10+ Reads ~ Assembled Reads (WAE Aln.)", xlab = "Assembled Reads", 
     ylab = "n loci with 10x depth or higher", pch = 1)
legend("bottomright", legend = c("Library1", "Library2"),
       col = c("black", "red"), 
       pch = 1)


##### ---------- Yellow Perch Alignment ---------- #####
getwd()
setwd("/Users/samjohnson/Documents/Sauger_102824/GeneticData/RadMapReport/RadMap_pflav/") 

# load in both tables
pflav <- read.table(file = "mappingreport_pflav_reads.txt", header = TRUE, sep = "\t", na.strings = "")
samp <- read.table(file = "sample_sampleLib.txt", header = TRUE, sep = "\t")

# merge by sample id, conserve all rows
radmap_pflav <- merge(samp, pflav, by = "sample", all = TRUE)

# the sampleLib column from samp contains both the plate and library info. also rename raw to rawreads and assembled to assembledreads
radmap_pflav <- radmap_pflav %>% 
  rename(plate = sampleLib) %>% 
  rename(rawreads = raw) %>% 
  rename(assembledreads = assembled)

# create a new column to JUST store the library, so that we can color by lib in downstream plots
radmap_pflav <- data.frame(radmap_pflav, lib = NA)

# rearrange columns so that plate and library info are at the left of the df
radmap_pflav <- radmap_pflav[, c("sample", "plate", "lib", "rawreads", "assembledreads", "lociN", "meanDepth", "lociNmin10reads", "meanDepthMin10reads")]

# for every row in the df, if the plate is from lib 1, put 1 in the lib column, else 2 for lib 2.
for(i in 1:nrow(radmap_pflav)) {
  if(grepl("Lib1", radmap_pflav$plate[i])) {
    radmap_pflav$lib[i] <- 1
  } else {
    radmap_pflav$lib[i] <- 2
  }
}

# Now the full dataframe has been established as we wish it to be.
setwd("/Users/samjohnson/Documents/Sauger_102824/GeneticData/RadMapReport/RadMap_pflav/CompleteReportPFlav/")
write.table(radmap_pflav, file = "mappingreport_pflav_reads_lib.txt")

plot(radmap_pflav$lociNmin10reads ~ radmap_pflav$assembledreads, col = radmap_pflav$lib,
     main = "Loci with 10+ Reads ~ Assembled Reads (YPE Aln.)", xlab = "Assembled Reads", 
     ylab = "n loci with 10x depth or higher", pch = 1)
legend("bottomright", legend = c("Library1", "Library2"),
       col = c("black", "red"), 
       pch = 1)

##### ---------- Plot Both Alignments Together  ---------- #####
#assign colors to libraries for each sp alignment
svit_col <- c("1" = "dodgerblue4", "2" = "deepskyblue1")
pflav_col <- c("1" = "tomato4", "2" = "salmon")

# create a color vector based on the 'lib' column for each sp alignment
svit_colors <- svit_col[radmap_svit$lib]
pflav_colors <- pflav_col[radmap_pflav$lib]

#now we'll use plot and points...
plot(radmap_svit$lociNmin10reads ~ radmap_svit$assembledreads, col = svit_colors,
     main = "Loci with 10+ Reads ~ Assembled Reads (by alignment)", xlab = "Assembled Reads", 
     ylab = "n loci with 10x depth or higher", pch = 19)
points(radmap_pflav$lociNmin10reads ~ radmap_pflav$assembledreads, col = radmap_pflav$lib, pch = 19)
legend("bottomright", legend = c("WAE_Lib1", "WAE_Lib2", "YPE_Lib1", "YPE_Lib2"),
       col = c("dodgerblue4", "deepskyblue1", "tomato4", "salmon"), 
       pch = 19)
