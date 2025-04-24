##################################################
## sauger_cross_check_trios.R by JPJ 3 vii 24
# modified by SPJ 042425
##################################################

## USAGE: Rscript sauger_cross_check_trios.R
## PURPOSE: to count proportion of correct trios inferred by sequoia
library(dplyr)
library(tidyverse)

nloci <- c(100, 500, 1000, 5000, 10000)
wild_samp <- c(0.2, 0.4, 0.6, 0.8, 1.0)
md <- c(0.01, 0.05, 0.1, 0.25, 0.5)

# Create parameter grid
param_grid <- expand.grid(nloci = nloci, wild_samp = wild_samp, md = md,
                          stringsAsFactors = FALSE)
out_trios <- param_grid %>%
  mutate(TP = 0, FP = 0, n_trios = 0)

for (i in seq_len(nrow(out_trios))) {
  nl <- out_trios$nloci[i]
  ws <- out_trios$wild_samp[i]
  miss <- out_trios$md[i]
  
  trio_file <- paste0("trios_", nl, "_", ws, "_", miss, ".txt")
  parent_file <- paste0("true_parents_", nl, "_", ws, "_", miss, ".txt")
		
  if (file.exists(trio_file) && file.exists(parent_file) &&
      file.info(trio_file)$size > 1 && file.info(parent_file)$size > 1) {
    
      trios <- read.table(trio_file, header = FALSE, sep = " ")
      true_parents <- read.table(parent_file, header = FALSE, sep = " ")
      
      # Loop over each inferred trio
      for (j in 1:nrow(trios)) {
        offspring <- trios[j, 1]
        inferred_parents <- trios[j, 2:3]
      
      # Get the true parents for this offspring
        true_pair <- true_parents[true_parents[, 1] == offspring, 2:3]
      
      # Compare sets (unordered) — full credit only if both parents match
        if (length(true_pair) == 2 &&
            setequal(as.character(inferred_parents), as.character(true_pair))) {
          out_trios$TP[i] <- out_trios$TP[i] + 1
        } else {
          out_trios$FP[i] <- out_trios$FP[i] + 1
        }
        
        out_trios$n_trios[i] <- out_trios$n_trios[i] + 1
      }
        
  } else {
    message(paste("Skipping due to missing or empty file:", trio_file, "or", parent_file))
    
    # Mark skipped cases with NA
    out_trios$TP[i] <- NA
    out_trios$FP[i] <- NA
    out_trios$n_trios[i] <- NA
  }
}
        
write.csv(out_trios, file = "tpfp_trios_md_sj.txt", sep = ",", row.names = FALSE)
        
# note, I just found some instances where F0 fish have F1 parents! Based solely on % DNA shared. 
# We need to try to specify generations I think...      


### --- Investigate Output --- ###
# NOTE: must comment out if we are to run the above portion again!!
# (first scp the outputfile to a local directory)

setwd("/Users/samjohnson/Desktop/")
trios <- read.csv(file = "tpfp_trios_md_sj.txt", header = TRUE, sep = ",")
trios <- data.frame(trios, truedivtot = NA, falsedivtot = NA, trueprop = NA, trueminusfalse = NA)
for (i in 1:nrow(trios)) {
  t <- trios$TP[i]
  f <- trios$FP[i]
  tot <- trios$n_trios[i]  # Was 'n_assignments' in your original code, but your file probably uses 'n_trios'
  
  if (!is.na(tot) && tot > 0) {
    trios$truedivtot[i] <- t / tot
    trios$falsedivtot[i] <- f / tot
  } else {
    trios$truedivtot[i] <- NA
    trios$falsedivtot[i] <- NA
  }
  
  trios$trueprop[i] <- if (!is.na(t)) t / 2000 else NA
  trios$trueminusfalse[i] <- if (!is.na(t) && !is.na(f)) t - f else NA
}

wild_samp_vals <- sort(unique(trios$wild_samp))
color_map <- c("purple", "blue", "forestgreen", "orange", "red")
names(color_map) <- wild_samp_vals

md_vals <- sort(unique(trios$md))
shape_map <- c(21, 22, 23, 24, 25)
names(shape_map) <- md_vals

nloci_vals <- sort(unique(trios$nloci))
color_map_loci <- c("firebrick4", "firebrick", "firebrick3", "firebrick2", "firebrick1")
names(color_map_loci) <- nloci_vals

trios$plot_col <- color_map[as.character(trios$wild_samp)]
trios$plot_pch <- shape_map[as.character(trios$md)]
trios$nloci_col <- color_map_loci[as.character(trios$nloci)]

par(mfrow = c(2, 2))

plot(trios$nloci, trios$truedivtot, col = trios$plot_col, pch = trios$plot_pch,
     log = "x", xlab = "nloci", ylab = "ntrue/nassign",
     main = "Prop. Assignments Made that are True",
     cex = 1.5)

plot(trios$nloci, trios$trueprop, col = trios$plot_col, pch = trios$plot_pch,
     log = "x", xlab = "nloci", ylab = "ntrue/2000",
     main = "Prop. of 2000 Possible True Assignments Made",
     cex = 1.5)

plot(trios$nloci, trios$trueminusfalse, col = trios$plot_col, pch = trios$plot_pch,
     log = "x", xlab = "nloci", ylab = "ntrue - nfalse",
     main = "n True Assignments - n False Assignments",
     ylim = c(-4000, 3000),
     cex = 1.5)
abline(h = 2000, col = "black", lty = 2, lwd = 2)

plot(trios$nloci, trios$trueminusfalse, col = trios$plot_col, pch = trios$plot_pch,
     log = "x", xlab = "nloci", ylab = "ntrue - nfalse",
     main = "n True Assignments - n False Assignments",
     ylim = c(0, 2200),
     cex = 1.5)
abline(h = 2000, col = "black", lty = 2, lwd = 2)

# REVISIT AT SOME POINT SOON TO MAKE PLOTS WITH WILD_SAMP AS X AXIS VAR
plot(trios$wild_samp, trios$TP, col = trios$nloci_col)












