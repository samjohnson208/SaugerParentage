##################################################
## sauger_cross_check_pars.R by JPJ 15 vii 24
# modified by SPJ 042125
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
        
        # Check if both files exist and are non-empty
        if (file.exists(true_file) && file.exists(pars_file) &&
            file.info(true_file)$size > 0 && file.info(pars_file)$size > 0) {
          
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
          message(paste("Skipping due to missing or empty file:", true_file, "or", pars_file))
        }
      }
    }
  }
}

# Rename columns for clarity and write output
colnames(out_pars) <- c("nloci", "wild_samp", "missing_data", "generation", "TP", "FP", "n_assignments")
write.csv(out_pars, file = "tpfp_pars_md_sj.txt", row.names = FALSE)


