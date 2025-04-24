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
        
write.csv(out_trios, file = "tftp_trios_md_sj.txt", sep = ",", row.names = FALSE)
        
# note, I just found some instances where F0 fish have F1 parents! Based solely on % DNA shared. 
# We need to try to specify generations I think...      
      
      
      
      