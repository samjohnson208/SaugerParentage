############################################
## sauger_cross_sim.R by JPJ 2 vii 24
############################################

## USAGE: Rscript sauger_cross_sim.R
## PURPOSE: Alternate simulation approach to slim

set.seed(17)


## set up starting parameters
args <- commandArgs(TRUE)             ## brings in values from stdin
nloci <- as.numeric(args[1])          ## number of loci
wild_samp <- as.numeric(args[2])      ## proportion of wild individuals (f0 non-hatchery parents, all f1s and f2s) sampled
miss <- as.numeric(args[3])           ## % of data that is missing per site

## fixed parameters
f0_inds <- 1000       ## total number of inds in F0 matrix
f0_cross_prop <- 0.2  ## proportion of f0 inds used to make hatchery crosses
f1h_inds <- 500      ## total number of inds in F1 hatchery matrix
f1w_inds <- 500      ## total number of inds in F1 wild matrix
f2_inds <- 1000       ## total number of inds in F2 matrix


## set up matrix of true parent ids
all_parents_full <- matrix(NA, (f0_inds+f1h_inds+f1w_inds+f2_inds), 3)


## set up f0 genotype matrix (NOTE!!!: this matrix is transposed relative to the others) and add to true parent ids matrix
f0_mat <- matrix(NA, nloci, f0_inds)
for (i in 1:nloci) {
	freq <- rbeta(1, 0.48, 0.206)  ## beta estimates come from fit to empirical data
	for (j in 1:f0_inds) {
		f0_mat[i,j] <- sum(rbinom(2, 1, freq))
		all_parents_full[j,1] <- paste0("f0_", j)
	}
}


## set up f1 hatchery matrix and add to true parent ids matrix
cross_inds <- floor(f0_cross_prop * f0_inds)
f1h_parents <- matrix(NA, f1h_inds, 2)
for (k in 1:f1h_inds) {
	f1h_parents[k,] <- sample(1:cross_inds, 2, replace=FALSE)
	all_parents_full[(k+f0_inds),1] <- paste0("f1h_", k)
	all_parents_full[(k+f0_inds),2] <- paste0("f0_", f1h_parents[k,1])
	all_parents_full[(k+f0_inds),3] <- paste0("f0_", f1h_parents[k,2])
}

f1h_mat <- matrix(NA, f1h_inds, nloci)
for (l in 1:f1h_inds) {
	for (m in 1:nloci) {
		if      (f0_mat[m,f1h_parents[l,1]]==0 && f0_mat[m,f1h_parents[l,2]]==0) { f1h_mat[l,m] <- 0 }
		else if (f0_mat[m,f1h_parents[l,1]]==0 && f0_mat[m,f1h_parents[l,2]]==1) { f1h_mat[l,m] <- rbinom(1, 1, 0.5) }
		else if (f0_mat[m,f1h_parents[l,1]]==0 && f0_mat[m,f1h_parents[l,2]]==2) { f1h_mat[l,m] <- 1 }
		else if (f0_mat[m,f1h_parents[l,1]]==1 && f0_mat[m,f1h_parents[l,2]]==0) { f1h_mat[l,m] <- rbinom(1, 1, 0.5) }
		else if (f0_mat[m,f1h_parents[l,1]]==1 && f0_mat[m,f1h_parents[l,2]]==1) { f1h_mat[l,m] <- sum(rbinom(2, 1, 0.5)) }
		else if (f0_mat[m,f1h_parents[l,1]]==1 && f0_mat[m,f1h_parents[l,2]]==2) { f1h_mat[l,m] <- rbinom(1, 1, 0.5) + 1 }
		else if (f0_mat[m,f1h_parents[l,1]]==2 && f0_mat[m,f1h_parents[l,2]]==0) { f1h_mat[l,m] <- 1 }
		else if (f0_mat[m,f1h_parents[l,1]]==2 && f0_mat[m,f1h_parents[l,2]]==1) { f1h_mat[l,m] <- rbinom(1, 1, 0.5) + 1 }
		else if (f0_mat[m,f1h_parents[l,1]]==2 && f0_mat[m,f1h_parents[l,2]]==2) { f1h_mat[l,m] <- 2 }
	}
}


## set up f1 wild matrix and add to true parent ids matrix
f1w_parents <- matrix(NA, f1w_inds, 2)
for (n in 1:f1w_inds) {
	f1w_parents[n,] <- sample((cross_inds+1):f0_inds, 2, replace=FALSE)
	all_parents_full[(n+f0_inds+f1h_inds),1] <- paste0("f1w_", n)
	all_parents_full[(n+f0_inds+f1h_inds),2] <- paste0("f0_", f1w_parents[n,1])
	all_parents_full[(n+f0_inds+f1h_inds),3] <- paste0("f0_", f1w_parents[n,2])
}

f1w_mat <- matrix(NA, f1w_inds, nloci)
for (o in 1:f1w_inds) {
	for (p in 1:nloci) {
		if      (f0_mat[p,f1w_parents[o,1]]==0 && f0_mat[p,f1w_parents[o,2]]==0) { f1w_mat[o,p] <- 0 }
		else if (f0_mat[p,f1w_parents[o,1]]==0 && f0_mat[p,f1w_parents[o,2]]==1) { f1w_mat[o,p] <- rbinom(1, 1, 0.5) }
		else if (f0_mat[p,f1w_parents[o,1]]==0 && f0_mat[p,f1w_parents[o,2]]==2) { f1w_mat[o,p] <- 1 }
		else if (f0_mat[p,f1w_parents[o,1]]==1 && f0_mat[p,f1w_parents[o,2]]==0) { f1w_mat[o,p] <- rbinom(1, 1, 0.5) }
		else if (f0_mat[p,f1w_parents[o,1]]==1 && f0_mat[p,f1w_parents[o,2]]==1) { f1w_mat[o,p] <- sum(rbinom(2, 1, 0.5)) }
		else if (f0_mat[p,f1w_parents[o,1]]==1 && f0_mat[p,f1w_parents[o,2]]==2) { f1w_mat[o,p] <- rbinom(1, 1, 0.5) + 1 }
		else if (f0_mat[p,f1w_parents[o,1]]==2 && f0_mat[p,f1w_parents[o,2]]==0) { f1w_mat[o,p] <- 1 }
		else if (f0_mat[p,f1w_parents[o,1]]==2 && f0_mat[p,f1w_parents[o,2]]==1) { f1w_mat[o,p] <- rbinom(1, 1, 0.5) + 1 }
		else if (f0_mat[p,f1w_parents[o,1]]==2 && f0_mat[p,f1w_parents[o,2]]==2) { f1w_mat[o,p] <- 2 }
	}
}


## merge the f1 matrices
f1_joined <- rbind(f1h_mat, f1w_mat)


## set up the f2 matrix and add to true parent ids matrix
f2_parents <- matrix(NA, f2_inds, 2)
for (q in 1:f2_inds) {
	f2_parents[q,] <- sample(1:dim(f1_joined)[1], 2, replace=FALSE)
	all_parents_full[(q+f0_inds+f1h_inds+f1w_inds),1] <- paste0("f2_", q)
	if (f2_parents[q,1] > f1h_inds) { all_parents_full[(q+f0_inds+f1h_inds+f1w_inds),2] <- paste0("f1w_", f2_parents[q,1]-f1h_inds) }
	else                            { all_parents_full[(q+f0_inds+f1h_inds+f1w_inds),2] <- paste0("f1h_", f2_parents[q,1]) }
	if (f2_parents[q,2] > f1h_inds) { all_parents_full[(q+f0_inds+f1h_inds+f1w_inds),3] <- paste0("f1w_", f2_parents[q,2]-f1h_inds) }
	else                            { all_parents_full[(q+f0_inds+f1h_inds+f1w_inds),3] <- paste0("f1h_", f2_parents[q,2]) }
}

f2_mat <- matrix(NA, f2_inds, nloci)
for (r in 1:f2_inds) {
	for (s in 1:nloci) {
		if      (f1_joined[f2_parents[r,1],s]==0 && f1_joined[f2_parents[r,2],s]==0) { f2_mat[r,s] <- 0 }
		else if (f1_joined[f2_parents[r,1],s]==0 && f1_joined[f2_parents[r,2],s]==1) { f2_mat[r,s] <- rbinom(1, 1, 0.5) }
		else if (f1_joined[f2_parents[r,1],s]==0 && f1_joined[f2_parents[r,2],s]==2) { f2_mat[r,s] <- 1 }
		else if (f1_joined[f2_parents[r,1],s]==1 && f1_joined[f2_parents[r,2],s]==0) { f2_mat[r,s] <- rbinom(1, 1, 0.5) }
		else if (f1_joined[f2_parents[r,1],s]==1 && f1_joined[f2_parents[r,2],s]==1) { f2_mat[r,s] <- sum(rbinom(2, 1, 0.5)) }
		else if (f1_joined[f2_parents[r,1],s]==1 && f1_joined[f2_parents[r,2],s]==2) { f2_mat[r,s] <- rbinom(1, 1, 0.5) + 1 }
		else if (f1_joined[f2_parents[r,1],s]==2 && f1_joined[f2_parents[r,2],s]==0) { f2_mat[r,s] <- 1 }
		else if (f1_joined[f2_parents[r,1],s]==2 && f1_joined[f2_parents[r,2],s]==1) { f2_mat[r,s] <- rbinom(1, 1, 0.5) + 1 }
		else if (f1_joined[f2_parents[r,1],s]==2 && f1_joined[f2_parents[r,2],s]==2) { f2_mat[r,s] <- 2 }
	}
}


## join all individuals together 
all_inds <- rbind(t(f0_mat), f1_joined, f2_mat)
all_inds_names <- matrix(NA, dim(all_inds)[1], 1)
for (t in 1:dim(all_inds)[1]) {
	if (t <= f0_inds) { all_inds_names[t,1] <- paste0("f0_", t) }
	else if (t <= f0_inds+f1h_inds) { all_inds_names[t,1] <- paste0("f1h_", t-f0_inds) }
	else if (t <= f0_inds+f1h_inds+f1w_inds) { all_inds_names[t,1] <- paste0("f1w_", t-f0_inds-f1h_inds) }
	else { all_inds_names[t,1] <- paste0("f2_", t-f0_inds-f1h_inds-f1w_inds) }
}
rownames(all_inds) <- all_inds_names


## remove a proportion of wild individuals (f0 non-hatchery parents, all f1s and f2s) based on wild_samp

if (wild_samp!=1) {
	wild_to_remove <- sample((cross_inds+1):dim(all_inds)[1], length((cross_inds+1):dim(all_inds)[1])-(length((cross_inds+1):dim(all_inds)[1])*wild_samp), replace=FALSE)
	final_inds <- all_inds[-c(wild_to_remove),]
	final_parents <- all_parents_full[-c(wild_to_remove),]
} else if (wild_samp==1) {
	final_inds <- all_inds
	final_parents <- all_parents_full
}

## remove a sample of x% from each site where x = md
for (u in 1:ncol(final_inds)) {
  n <- nrow(final_inds)
  sample_indices <- sample(n, size = floor(miss * n), replace = FALSE)
  final_inds[sample_indices, u] <- -9
}

## sequoia
library(sequoia)
output <- GetMaybeRel(final_inds)


## write out files
write.table(output$MaybeTrio, file=paste0("trios_", nloci, "_", wild_samp, "_", miss, ".txt"), row.names=F, quote=F)
write.table(output$MaybePar, file=paste0("pars_", nloci, "_", wild_samp, "_", miss, ".txt"), row.names=F, quote=F)
write.table(final_parents, file=paste0("true_parents_", nloci, "_", wild_samp, "_", miss, ".txt"), row.names=F, col.names=F, quote=F)