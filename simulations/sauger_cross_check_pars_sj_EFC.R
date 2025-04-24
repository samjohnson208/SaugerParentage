##################################################
## sauger_cross_check_pars.R by JPJ 15 vii 24
# modified by SPJ 042125 to sauger_cross_check_pars_sj_EFC.R 
# where EFC stands for empty file check!
##################################################

## USAGE: Rscript sauger_cross_check_pars.R
## PURPOSE: to count proportion of correct parent-offspring relationships inferred by sequoia

library(tidyverse)
library(dplyr)

# Define simulation parameter vectors
nloci <- c(100, 500, 1000, 5000, 10000)
wild_samp <- c(0.2, 0.4, 0.6, 0.8, 1.0)
md <- c(0.01, 0.05, 0.1, 0.25, 0.5)
gen_cats <- c("f0", "f1h", "f1w", "f2")

# Create output data frame to store results
param_grid <- expand.grid(nloci = nloci, wild_samp = wild_samp, md = md, gen_cat = gen_cats,
                          stringsAsFactors = FALSE)
out_pars <- param_grid %>%
  mutate(TP = 0, FP = 0, n_assignments = 0)

# Main loop over all parameter combinations
for (g in 1:length(nloci)) {
  for (h in 1:length(wild_samp)) {
    for (z in 1:length(md)) {
      for (f in 1:length(gen_cats)) {
        
        gen <- gen_cats[f]
        nl <- nloci[g]
        ws <- wild_samp[h]
        miss <- md[z]
        
        true_file <- paste0("true_parents_", nl, "_", ws, "_", miss, ".txt")
        pars_file <- paste0("pars_", nl, "_", ws, "_", miss, ".txt")
        
        # Check if both files exist and are non-empty (not 1 byte)
        if (file.exists(true_file) && file.exists(pars_file) &&
            file.info(true_file)$size > 1 && file.info(pars_file)$size > 1) {
          
          truth <- read.table(true_file, header = FALSE, sep = " ")
          pars <- read.table(pars_file, header = TRUE, sep = " ")
          
          for (i in 1:nrow(pars)) {
            offspring_id <- pars[i, 2]
            parent_id <- pars[i, 1]
            
            row_index <- which(out_pars$nloci == nl &
                                 out_pars$wild_samp == ws &
                                 out_pars$md == miss &
                                 out_pars$gen_cat == gen)
            
            out_pars$n_assignments[row_index] <- out_pars$n_assignments[row_index] + 1
            
            if (grepl("f0", offspring_id)) {
              out_pars$FP[row_index] <- out_pars$FP[row_index] + 1
            } else {
              true_parents <- truth[truth[, 1] == offspring_id, 2:3]
              if (parent_id %in% as.character(unlist(true_parents))) {
                out_pars$TP[row_index] <- out_pars$TP[row_index] + 1
              } else {
                out_pars$FP[row_index] <- out_pars$FP[row_index] + 1
              }
            }
          }
        } else {
          # Skip empty or nearly empty files and set NA for results
          message(paste("Skipping due to missing or empty file:", true_file, "or", pars_file))
          
          row_index <- which(out_pars$nloci == nl &
                               out_pars$wild_samp == ws &
                               out_pars$md == miss &
                               out_pars$gen_cat == gen)
          
          # Assign NA to TP, FP, and n_assignments for these combinations
          out_pars$TP[row_index] <- NA
          out_pars$FP[row_index] <- NA
          out_pars$n_assignments[row_index] <- NA
        }
      }
    }
  }
}

# Rename columns for clarity and write output
colnames(out_pars) <- c("nloci", "wild_samp", "missing_data", "generation", "TP", "FP", "n_assignments")
write.csv(out_pars, file = "tpfp_pars_md_sj.csv", row.names = FALSE)


### --- Investigate Output --- ###
# NOTE: must comment out if we are to run the above portion again!!
# (first scp the outputfile to a local directory)
setwd("/Users/samjohnson/Desktop/")
pars <- read.csv(file = "tpfp_pars_md_sj.txt", header = TRUE, sep = ",")
pars <- data.frame(pars, truedivtot = NA, falsedivtot = NA, trueprop = NA, trueminusfalse = NA)
for(i in 1:nrow(pars)){
  t <- pars$TP[i]
  f <- pars$FP[i]
  tot <- pars$n_assignments[i]
  pars$truedivtot[i] <- t/tot
  pars$falsedivtot[i] <- f/tot
  pars$trueprop[i] <- t/4000
  pars$trueminusfalse[i] <- t-f
}

wild_samp_vals <- sort(unique(pars$wild_samp))
color_map <- c("purple", "blue", "forestgreen", "orange", "red")
names(color_map) <- wild_samp_vals

md_vals <- sort(unique(pars$missing_data))
shape_map <- c(21, 22, 23, 24, 25)
names(shape_map) <- md_vals

nloci_vals <- sort(unique(pars$nloci))
color_map_loci <- c("firebrick4", "firebrick", "firebrick3", "firebrick2", "firebrick1")
names(color_map_loci) <- nloci_vals

pars$plot_col <- color_map[as.character(pars$wild_samp)]
pars$plot_pch <- shape_map[as.character(pars$missing_data)]
pars$nloci_col <- color_map_loci[as.character(pars$nloci)]

par(mfrow = c(2, 2))

plot(pars$nloci, pars$truedivtot, col = pars$plot_col, pch = pars$plot_pch,
     log = "x", xlab = "nloci", ylab = "ntrue/nassign",
     main = "Prop. Assignments Made that are True",
     cex = 1.5)

plot(pars$nloci, pars$trueprop, col = pars$plot_col, pch = pars$plot_pch,
     log = "x", xlab = "nloci", ylab = "ntrue/4000",
     main = "Prop. of 4000 Possible True Assignments Made",
     cex = 1.5)

plot(pars$nloci, pars$trueminusfalse, col = pars$plot_col, pch = pars$plot_pch,
     log = "x", xlab = "nloci", ylab = "ntrue - nfalse",
     main = "n True Assignments - n False Assignments",
     ylim = c(-20000, 5000),
     cex = 1.5)
abline(h = 4000, col = "black", lty = 2, lwd = 2)

plot(pars$nloci, pars$trueminusfalse, col = pars$plot_col, pch = pars$plot_pch,
     log = "x", xlab = "nloci", ylab = "ntrue - nfalse",
     main = "n True Assignments - n False Assignments",
     ylim = c(0, 5000),
     cex = 1.5)
abline(h = 4000, col = "black", lty = 2, lwd = 2)

# REVISIT AT SOME POINT SOON TO MAKE PLOTS WITH WILD_SAMP AS X AXIS VAR
plot(pars$wild_samp, pars$TP, col = pars$nloci_col)
