##################################################
## sauger_cross_check_pars.R by JPJ 15 vii 24
# modified by SPJ 042125
##################################################

## USAGE: Rscript sauger_cross_check_pars.R
## PURPOSE: to count proportion of correct parent-offspring relationships inferred by sequoia
library(tidyverse)
library(dplyr)

nloci <- c(100, 500, 1000, 5000, 10000)
wild_samp <- c(0.2, 0.4, 0.6, 0.8, 1.0)
md <- c(0.01, 0.05, 0.1, 0.25, 0.5)

gen_cats <- c("f0", "f1h", "f1w", "f2")

## set up output matrix - JPJ
# out_pars <- matrix(0, length(nloci)*length(wild_samp)*length(md)*length(gen_cats), 5)
# colnames(out_pars) <- c("nloci", "wild_samp", "md", "gen_cat", "TP", "FP")


### --- SPJ Modification of Output Matrix and Looping Structure --- ###
param_grid <- expand.grid(nloci = nloci, wild_samp = wild_samp, md = md, gen_cat = gen_cats,
                          stringsAsFactors = FALSE)
out_pars <- param_grid %>%
  mutate(TP = 0, FP = 0, n_assign = 0)

ctr <- 0
for (d in 1:length(nloci)) {
  for (e in 1:length(wild_samp)) {
    for (c in 1:length(md)) {
      for (f in 1:length(gen_cats)) {
        ctr <- ctr + 1
        out_pars[ctr, 1] <- nloci[d]
        out_pars[ctr, 2] <- wild_samp[e]
        out_pars[ctr, 3] <- md[c]
        out_pars[ctr, 4] <- gen_cats[f]
        out_pars[ctr, 5] <- 0  # TP
        out_pars[ctr, 6] <- 0  # FP
        out_pars[ctr, 7] <- 0  # n_assignments
      }
    }
  }
}

for (g in 1:length(nloci)) {
  for (h in 1:length(wild_samp)) {
    for (z in 1:length(md)) {
      
      truth_file <- paste0("true_parents_", nloci[g], "_", wild_samp[h], "_", md[z], ".txt")
      pars_file <- paste0("pars_", nloci[g], "_", wild_samp[h], "_", md[z], ".txt")
      
      truth <- read.delim(truth_file, header = FALSE, sep = " ")
      pars <- read.delim(pars_file, header = TRUE, sep = " ")
      
      for (i in 1:nrow(pars)) {
        offspring_id <- pars[i, 2]
        parent_id <- pars[i, 1]
        
        gen_cat <- if (grepl("f0", offspring_id)) {
          "f0"
        } else if (grepl("f1h", offspring_id)) {
          "f1h"
        } else if (grepl("f1w", offspring_id)) {
          "f1w"
        } else if (grepl("f2", offspring_id)) {
          "f2"
        } else {
          next
        }
        
        match_row <- out_pars[,1] == nloci[g] &
          out_pars[,2] == wild_samp[h] &
          out_pars[,3] == md[z] &
          out_pars[,4] == gen_cat
        
        true_parents <- truth[truth[,1] == offspring_id, 2:3]
        
        out_pars[match_row, 7] <- as.numeric(out_pars[match_row, 7]) + 1  # increment n_assignments
        
        if (nrow(true_parents) > 0 && parent_id %in% unlist(true_parents)) {
          out_pars[match_row, 5] <- as.numeric(out_pars[match_row, 5]) + 1  # TP
        } else {
          out_pars[match_row, 6] <- as.numeric(out_pars[match_row, 6]) + 1  # FP
        }
      }
    }
  }
}

# Save results
colnames(out_pars) <- c("nloci", "wild_samp", "missing_data", "generation", "TP", "FP", "n_assignments")
write.csv(out_pars, file = "tpfp_pars_md_sj.txt", row.names = FALSE)


# to do when we come back
# create a folder with only the files with md
# see how we ran this the first time (r script nested in a slurm?) i actually don't think i did that so 
# i'll have to make a slurm?
# ensure it works. 
# then do trios.


### --- JPJ Original Looping Portion --- ###

# ## read in sequoia parents output and count true and false positives (TP & FP)
# for (g in 1:length(nloci)) {
# 	for (h in 1:length(wild_samp)) {
# 	  for (z in 1:length(md)) {
# 		# truth <- read.delim(paste0("true_parents_", nloci[g], "_", wild_samp[h], ".txt"), header=FALSE, sep=" ") ## true parental relationships
# 		# pars <- read.delim(paste0("pars_", nloci[g], "_", wild_samp[h], ".txt"), header=TRUE, sep=" ")	 ## parental relationships from sequoia
# 		for (i in 1:dim(pars)[1]) {
# 			if (grepl("f0", pars[i,2])==TRUE) { ## if offspring is from f0 generation
# 				out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f0", 5] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f0", 5]) + 1
# 			}
# 			else if (grepl("f1h", pars[i,2])==TRUE) { ## if offspring is from f1h generation
# 				if (pars[i,1] %in% truth[truth[,1]==pars[i,2],2:3]==TRUE) { out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1h", 4] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1h", 4]) + 1 }
# 				else if (pars[i,1] %in% truth[truth[,1]==pars[i,2],2:3]==FALSE) { out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1h", 5] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1h", 5]) + 1 }
# 			}
# 			else if (grepl("f1w", pars[i,2])==TRUE) { ## if offspring is from f1w generation
# 				if (pars[i,1] %in% truth[truth[,1]==pars[i,2],2:3]==TRUE) { out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1w", 4] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1w", 4]) + 1 }
# 				else if (pars[i,1] %in% truth[truth[,1]==pars[i,2],2:3]==FALSE) { out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1w", 5] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1w", 5]) + 1 }
# 			}
# 			else if (grepl("f2", pars[i,2])==TRUE) { ## if offspring is from f2 generation
# 				if (pars[i,1] %in% truth[truth[,1]==pars[i,2],2:3]==TRUE) { out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f2", 4] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f2", 4]) + 1 }
# 				else if (pars[i,1] %in% truth[truth[,1]==pars[i,2],2:3]==FALSE) { out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f2", 5] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f2", 5]) + 1 }
# 
# 		  	}
# 	  	}
#   	}
# 	}
# }
# 
# #write.table(out_pars, file="tpfp_pars.txt", row.names=FALSE, quote=FALSE)