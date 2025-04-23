##################################################
## sauger_cross_check_trios.R by JPJ 3 vii 24
# modified by SPJ 042125
##################################################

## USAGE: Rscript sauger_cross_check_trios.R
## PURPOSE: to count proportion of correct trios inferred by sequoia

nloci <- c(100, 500, 1000, 5000, 10000)
wild_samp <- c(0.2, 0.4, 0.6, 0.8, 1.0)
md <- c(0.01, 0.05, 0.1, 0.25, 0.5)

out_trios <- matrix(NA, length(nloci), length(wild_samp))
rownames(out_trios) <- nloci
colnames(out_trios) <- wild_samp

for (g in 1:length(nloci)) {
	for (h in 1:length(wild_samp)) {
		## read in true parent ids
		true_parents <- read.delim(paste0("true_parents_", nloci[g], "_", wild_samp[h], ".txt"), header=FALSE, sep=" ")
		## read in inferred relationships by sequoia and check accuracy
		trios <- read.delim(paste0("trios_", nloci[g], "_", wild_samp[h], ".txt"), header=TRUE, sep=" ")
		right_count <- 0
		wrong_count <- 0
		for (i in 1:dim(trios)[1]) {
			if (trios[i,2] %in% true_parents[true_parents[,1]==trios[i,1],2:3]==TRUE) { right_count <- right_count+1 }
			else                                                                      { wrong_count <- wrong_count+1 }
			if (trios[i,3] %in% true_parents[true_parents[,1]==trios[i,1],2:3]==TRUE) { right_count <- right_count+1 }
			else                                                                      { wrong_count <- wrong_count+1 }
		}
		out_trios[g,h] <- right_count/(dim(trios)[1]*2)
	}
}

write.table(out_trios, file="true_parents_in_trios.txt", quote=FALSE)