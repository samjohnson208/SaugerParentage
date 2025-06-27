## hist_meandepth.R by SPJ 062625
## PURPOSE: generate a histogram of mean depth of aligned and initially filtered
# reads to see what min/max meandepth filter makes reasonable sense.
## USAGE: Rscript hist_meandepth.R

## GENERATING INPUT FILES RELIES ON SaugerParentage/slurms/slurm_sauger_site_mean_depth.sh ## 


setwd("/Users/samjohnson/Desktop/")
smd <- read.table(file = "variants_pflav_mem_t2_bial_noindels_q20.ldepth.mean", header = TRUE)
colnames(smd)
hist(smd$MEAN_DEPTH, breaks = 20, xlim = c(0,40))
summary(smd$MEAN_DEPTH)

