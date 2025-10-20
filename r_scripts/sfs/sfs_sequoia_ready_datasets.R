# sfs_sequoia_ready_datasets.R by SPJ 102025
# purpose: generate sfs for all 17 datasets that we have prepared for sequoia so
# far. see dataset prep in Sequoia_ContamFilt_mindep8_md5_RADseq.R

# i have these datasets loaded into my environment from that script, so i will
# attempt to loop over all of them in this script.

# all this comes from a dataset that has been...
# cleaned using illumina, phix, and ecoli databases, and has been cleaned using
# FASTP (remove bases with q<20, reads shorter than 50bp, poly g tails, and any 
# remaining sequence adapters). This dataset was then aligned to the walleye reference
# using BWA MEM, variant calling was run, vcf was reheadered, and the data were 
# filtered using the slurm_1 and slurm_3 filtering scripts for the following params.

# rawfiltered - see bcftools mpileup in slurm_sauger_variants.sh
# contam - contaminant filtering mentioned above
# fastp - see above
# bial - keep only biallelic sites
# no indels - keep sites that are not insertions or deletions
# q40 - keep sites with site quality > 40
# mindep8 - minimum mean read depth 4 or above (removes low quality sites, low confidence calls)
# maxdep75 - maximum mean read depth 75 or lower (removes paralogs)
# maf30 - minor allele frequency 30% or higher
# miss95 - each snp needs to have 5% or less missing data

# # and thinned to various degrees to filter for LD (25K, 50K, 75K, 100K to 900K, 
# 1M to 5M), for a total of 17 sets.

# each object is called, at this point: check_thin___K/M (these are input for GetMaybeRel())

# create a list of each of the objects
GenoM_list <- list(
  check_thin25K = check_thin25K,
  check_thin50K = check_thin50K,
  check_thin75K = check_thin75K,
  check_thin100K = check_thin100K,
  check_thin200K = check_thin200K,
  check_thin300K = check_thin300K,
  check_thin400K = check_thin400K,
  check_thin500K = check_thin500K,
  check_thin600K = check_thin600K,
  check_thin700K = check_thin700K,
  check_thin800K = check_thin800K,
  check_thin900K = check_thin900K,
  check_thin1M = check_thin1M,
  check_thin2M = check_thin2M,
  check_thin3M = check_thin3M,
  check_thin4M = check_thin4M,
  check_thin5M = check_thin5M
)

# set up plotting area
par(mfrow = c(5, 4), mar = c(4, 4, 2, 1))

# define shared 
n_breaks <- 50
xlims_unfold <- c(0, 1)
xlims_fold <- c(0, 0.5)

# -----------------------
# UNFOLDED SFS
# -----------------------

# set up plotting area
par(mfrow = c(5, 4), mar = c(4, 4, 2, 1))

for (name in names(GenoM_list)) {
  
  geno_mat <- GenoM_list[[name]]
  
  # compute site means (per SNP column)
  site_mean <- apply(geno_mat, 2, mean, na.rm = TRUE)
  
  # divide by two to get unfolded allele frequency
  site_mean_d2 <- site_mean / 2
  
  # plot unfolded SFS
  hist(
    site_mean_d2,
    breaks = n_breaks,
    main = paste("Unfolded SFS,", name),
    xlab = "Ref Allele Frequency",
    xlim = xlims_unfold,
    col = "royalblue",
    border = "white"
  )
}


# -----------------------
# FOLDED SFS
# -----------------------

# set up plotting area
par(mfrow = c(5, 4), mar = c(4, 4, 2, 1))

for (name in names(GenoM_list)) {
  
  geno_mat <- GenoM_list[[name]]
  
  # compute site means (per SNP column)
  site_mean <- apply(geno_mat, 2, mean, na.rm = TRUE)
  site_mean_d2 <- site_mean / 2
  
  # fold around 0.5
  fold_site_mean_d2 <- ifelse(site_mean_d2 > 0.5, 1 - site_mean_d2, site_mean_d2)
  
  # plot folded SFS
  hist(
    fold_site_mean_d2,
    breaks = n_breaks,
    main = paste("Folded SFS,", name),
    xlab = "Ref Allele Frequency",
    xlim = xlims_fold,
    col = "lightcoral",
    border = "white"
  )
}
