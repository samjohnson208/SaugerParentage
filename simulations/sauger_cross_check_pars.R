##################################################
## sauger_cross_check_pars.R by JPJ 15 vii 24
##################################################

## USAGE: Rscript sauger_cross_check_pars.R
## PURPOSE: to count proportion of correct parent-offspring relationships inferred by sequoia

nloci <- c(100, 500, 1000, 5000, 10000)
wild_samp <- c(0.2, 0.4, 0.6, 0.8, 1.0)
gen_cats <- c("f0", "f1h", "f1w", "f2")

## set up output matrix
out_pars <- matrix(0, length(nloci)*length(wild_samp)*length(gen_cats), 5)
colnames(out_pars) <- c("nloci", "wild_samp", "gen_cat", "TP", "FP")
ctr <- 0
for (d in 1:length(nloci)) {
	for (e in 1:length(wild_samp)) {
		for (f in 1:length(gen_cats)) {
			ctr <- ctr + 1
			out_pars[ctr,1] <- nloci[d]
			out_pars[ctr,2] <- wild_samp[e]
			out_pars[ctr,3] <- gen_cats[f]
		}
	}
}


## read in sequoia parents output and count true and false positives (TP & FP)
for (g in 1:length(nloci)) {
	for (h in 1:length(wild_samp)) {
		truth <- read.delim(paste0("true_parents_", nloci[g], "_", wild_samp[h], ".txt"), header=FALSE, sep=" ") ## true parental relationships
		pars <- read.delim(paste0("pars_", nloci[g], "_", wild_samp[h], ".txt"), header=TRUE, sep=" ")	 ## parental relationships from sequoia
		for (i in 1:dim(pars)[1]) {
			if (grepl("f0", pars[i,2])==TRUE) { ## if offspring is from f0 generation
				out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f0", 5] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f0", 5]) + 1
			}
			else if (grepl("f1h", pars[i,2])==TRUE) { ## if offspring is from f1h generation
				if (pars[i,1] %in% truth[truth[,1]==pars[i,2],2:3]==TRUE) { out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1h", 4] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1h", 4]) + 1 }
				else if (pars[i,1] %in% truth[truth[,1]==pars[i,2],2:3]==FALSE) { out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1h", 5] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1h", 5]) + 1 }
			}
			else if (grepl("f1w", pars[i,2])==TRUE) { ## if offspring is from f1w generation
				if (pars[i,1] %in% truth[truth[,1]==pars[i,2],2:3]==TRUE) { out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1w", 4] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1w", 4]) + 1 }
				else if (pars[i,1] %in% truth[truth[,1]==pars[i,2],2:3]==FALSE) { out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1w", 5] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f1w", 5]) + 1 }
			}
			else if (grepl("f2", pars[i,2])==TRUE) { ## if offspring is from f2 generation
				if (pars[i,1] %in% truth[truth[,1]==pars[i,2],2:3]==TRUE) { out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f2", 4] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f2", 4]) + 1 }
				else if (pars[i,1] %in% truth[truth[,1]==pars[i,2],2:3]==FALSE) { out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f2", 5] <- as.numeric(out_pars[out_pars[,1]==nloci[g] & out_pars[,2]==wild_samp[h] & out_pars[,3]=="f2", 5]) + 1 }

			}
		}
	}
}

write.table(out_pars, file="tpfp_pars.txt", row.names=FALSE, quote=FALSE)