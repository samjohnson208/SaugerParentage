# here's the goal. these are just notes okay? chill.

# we're keeping that same genom and errm and extrapolating the rest of the data.

# we want to run sequoia() on the pairwise gens, do the pedigree reconstruction
seq_f0_f2 <- sequoia(GenoM = check_thin100K_f0f2,
                     LifeHistData = LH_F0_F2,
                     args.AP=list(Discrete = TRUE, MaxAgeParent = 2),
                     Module = "ped",
                     Err = errM,
                     Complex = "full",
                     Herm = "no",
                     UseAge = "yes",
                     CalcLLR = TRUE,
                     StrictGenoCheck = TRUE,
                     DummyPrefix = c("F", "M"),
                     Tfilter = -2,
                     Tassign = 1.0)

# no validation necessary here from lookup tables? or do we want to do the validation
# for the f1 -> f0 assignments? 

# f0-f0, f0-f1, f1-f1
# f1-f2, f2-f2, f0-f2

# then we do another one where we do f0-f1, where we take the trios and do valid
# is true, false, and then all other f0-f1 toprels?

# okay, so for each one of those six things above, we should have dists for those
# relationships arranged in two panels. 1 should be all combinations, 2 should be
# the pairs that get assigned by sequoia()

setwd("/Users/samjohnson/Documents/Sauger_102325/GeneticData/Sequoia/Sequoia_Inp/contam_fastp_svit_mem/firstfilt_hardfilt_thinned/mindep8_maf30/geno_mat")
getwd()
inds_all <- read.csv(file = "Inds_F0_F1Spawn_F1Juv_F2.csv", header = TRUE)
table(inds_all$Group)

##### ---- run sequoia() for all pairs of gens ---- #####

# filter genotype matrix for all six runs
  # should have the code for this already. find those lh dfs.
dim(check_thin100K_all) #1030 inds. 943 snps.
# inds_all has 1077 inds. (everyone except for test f1's)
dim(LH_All) # also 1030 inds. where are those other 47 coming from?
# 17 filtered out cuz they didn't get sequenced, 
  # inds_all was build from SAR_Data, includes empty tubes that got id's
#30 filtered out for md per ind
  # (see Sequoia_ContamFilt_mindep8_md5_RADseqErr_ALLINDS.R)
inds_all <- inds_all %>% 
    filter(ID %in% rownames(check_thin100K_all))
dim(inds_all) # 1030 nice.
table(inds_all$Group)

##### ---- make filtered genotype matrices ---- #####
# f0-f0
f0_inds <- inds_all %>% 
    filter(Group == "F0")
check_thin100K_f0 <- check_thin100K_all[rownames(check_thin100K_all) %in% f0_inds$ID, , drop = FALSE]
dim(check_thin100K_f0) # 250 inds, 943 snps

# f1-f1
f1_inds <- inds_all %>% 
    filter(Group %in% c("F1 - Juvenile", "F1 - Spawning"))
check_thin100K_f1 <- check_thin100K_all[rownames(check_thin100K_all) %in% f1_inds$ID, , drop = FALSE]
dim(check_thin100K_f1) # 326 inds, 943 snps

# f2-f2
f2_inds <- inds_all %>% 
  filter(Group == "F2")
check_thin100K_f2 <- check_thin100K_all[rownames(check_thin100K_all) %in% f2_inds$ID, , drop = FALSE]
dim(check_thin100K_f2) # 454 inds, 943 snps

      # 250 + 326 + 454 = 1030 YEP.

# f0-f1
f0_f1_inds <- rbind(f0_inds, f1_inds)
dim(f0_f1_inds) # 576 inds
check_thin100K_f0f1 <- rbind(check_thin100K_f0, check_thin100K_f1)
dim(check_thin100K_f0f1) # 576 inds, 943 snps

# f1-f2 
f1_f2_inds <- rbind(f1_inds, f2_inds)
dim(f1_f2_inds) # 780 inds
check_thin100K_f1f2 <- rbind(check_thin100K_f1, check_thin100K_f2)
dim(check_thin100K_f1f2) # 780 inds, 943 snps

# f0-f2
f0_f2_inds <- rbind(f0_inds, f2_inds)
dim(f0_f2_inds) # 704 inds
check_thin100K_f0f2 <- rbind(check_thin100K_f0, check_thin100K_f2)
dim(check_thin100K_f0f2) # 704 inds, 943 snps

##### ---- make filtered LH dataframes ---- #####
  # should have the code for this already. find those lh dfs.
dim(LH_All) # 1030 inds.

LH_f0 <- LH_All %>% 
    filter(ID %in% rownames(check_thin100K_f0))
dim(LH_f0) # 250 inds

LH_f1 <- LH_All %>% 
  filter(ID %in% rownames(check_thin100K_f1))
dim(LH_f1) # 326 inds

LH_f2 <- LH_All %>% 
  filter(ID %in% rownames(check_thin100K_f2))
dim(LH_f2) # 454 inds

LH_f0f1 <- rbind(LH_f0, LH_f1)
dim(LH_f0f1) # 576 inds

LH_f1f2 <- rbind(LH_f1, LH_f2)
dim(LH_f1f2) # 780 inds

LH_f0f2 <- rbind(LH_f0, LH_f2)
dim(LH_f0f2) # 704 inds

##### ---- run sequoia(), GetMaybeRel(), GetRelM(), PlotRelPairs() ---- #####
  # do we run it for f0-f0 with the ageprior? can it detect all of the FS
  # if it's also searching for f0-f1 relationships? 

errM <- Err_RADseq(E0 = 0.075, E1 = 0.025, Return = 'matrix')

##### ---- f0_f0 ---- #####
seq_f0 <- sequoia(GenoM = check_thin100K_f0,
                  LifeHistData = LH_f0,
                  Module = "ped",
                  Err = errM,
                  Complex = "full",
                  Herm = "no",
                  UseAge = "yes",
                  args.AP=list(Discrete = TRUE, 
                               MinAgeParent = 1, MaxAgeParent = 1),
                  CalcLLR = TRUE,
                  StrictGenoCheck = TRUE,
                  DummyPrefix = c("F", "M"),
                  Tfilter = -2,
                  Tassign = 1.0)

gmr_f0 <- GetMaybeRel(GenoM = check_thin100K_f0,
                      AgePrior = seq_f0[["AgePriors"]],
                      Err = errM,
                      Module = "ped",
                      Complex = "full",
                      LifeHistData = LH_f0,
                      Herm = "no",
                      quiet = FALSE,
                      Tfilter = -2,
                      Tassign = 1.0,
                      MaxPairs = 7 * nrow(check_thin100K_f0))

relm_f0 <- GetRelM(Pedigree = seq_f0[["Pedigree"]],
                   GenBack = 1, 
                   patmat = FALSE,
                   directed = TRUE,
                   Return = 'Matrix')
table(unique(relm_f0))

relmf0_plot <- PlotRelPairs(RelM = relm_f0, pch.symbols = TRUE)

seq_f0_summary <- SummarySeq(SeqList = seq_f0)

##### ---- f0_f1 ---- #####
seq_f0f1 <- sequoia(GenoM = check_thin100K_f0f1,
                    LifeHistData = LH_f0f1,
                    Module = "ped",
                    Err = errM,
                    Complex = "full",
                    Herm = "no",
                    UseAge = "yes",
                    args.AP=list(Discrete = TRUE, 
                                 MinAgeParent = 1, MaxAgeParent = 1),
                    CalcLLR = TRUE,
                    StrictGenoCheck = TRUE,
                    DummyPrefix = c("F", "M"),
                    Tfilter = -2,
                    Tassign = 1.0)

gmr_f0f1 <- GetMaybeRel(GenoM = check_thin100K_f0f1,
                      AgePrior = seq_f0f1[["AgePriors"]],
                      Err = errM,
                      Module = "ped",
                      Complex = "full",
                      LifeHistData = LH_f0f1,
                      Herm = "no",
                      quiet = FALSE,
                      Tfilter = -2,
                      Tassign = 1.0,
                      MaxPairs = 7 * nrow(check_thin100K_f0f1))

relm_f0f1 <- GetRelM(Pedigree = seq_f0f1[["Pedigree"]],
                   GenBack = 1, 
                   patmat = FALSE,
                   directed = TRUE,
                   Return = 'Matrix')
table(unique(relm_f0f1))

relmf0f1_plot <- PlotRelPairs(RelM = relm_f0f1, drop.U = TRUE, pch.symbols = FALSE)

seq_f0f1_summary <- SummarySeq(SeqList = seq_f0f1)

##### ---- f1_f1 ---- #####

seq_f1 <- sequoia(GenoM = check_thin100K_f1,
                  LifeHistData = LH_f1,
                  Module = "ped",
                  Err = errM,
                  Complex = "full",
                  Herm = "no",
                  UseAge = "yes",
                  args.AP=list(Discrete = TRUE, 
                               MinAgeParent = 1, MaxAgeParent = 1),
                  CalcLLR = TRUE,
                  StrictGenoCheck = TRUE,
                  DummyPrefix = c("F", "M"),
                  Tfilter = -2,
                  Tassign = 1.0)

gmr_f1 <- GetMaybeRel(GenoM = check_thin100K_f1,
                      AgePrior = seq_f1[["AgePriors"]],
                      Err = errM,
                      Module = "ped",
                      Complex = "full",
                      LifeHistData = LH_f1,
                      Herm = "no",
                      quiet = FALSE,
                      Tfilter = -2,
                      Tassign = 1.0,
                      MaxPairs = 7 * nrow(check_thin100K_f1))

relm_f1 <- GetRelM(Pedigree = seq_f1[["Pedigree"]],
                   GenBack = 1, 
                   patmat = FALSE,
                   directed = TRUE,
                   Return = 'Matrix')
table(unique(relm_f1))

relmf1_plot <- PlotRelPairs(RelM = relm_f1, pch.symbols = TRUE)

seq_f1_summary <- SummarySeq(SeqList = seq_f1)


##### ---- f1_f2 ---- #####
seq_f1f2 <- sequoia(GenoM = check_thin100K_f1f2,
                    LifeHistData = LH_f1f2,
                    Module = "ped",
                    Err = errM,
                    Complex = "full",
                    Herm = "no",
                    UseAge = "yes",
                    args.AP=list(Discrete = TRUE, 
                                 MinAgeParent = 1, MaxAgeParent = 1),
                    CalcLLR = TRUE,
                    StrictGenoCheck = TRUE,
                    DummyPrefix = c("F", "M"),
                    Tfilter = -2,
                    Tassign = 1.0)

gmr_f1f2 <- GetMaybeRel(GenoM = check_thin100K_f1f2,
                        AgePrior = seq_f1f2[["AgePriors"]],
                        Err = errM,
                        Module = "ped",
                        Complex = "full",
                        LifeHistData = LH_f1f2,
                        Herm = "no",
                        quiet = FALSE,
                        Tfilter = -2,
                        Tassign = 1.0,
                        MaxPairs = 7 * nrow(check_thin100K_f1f2))

relm_f1f2 <- GetRelM(Pedigree = seq_f1f2[["Pedigree"]],
                     GenBack = 1, 
                     patmat = FALSE,
                     directed = TRUE,
                     Return = 'Matrix')
table(unique(relm_f1f2))

relmf1f2_plot <- PlotRelPairs(RelM = relm_f1f2, drop.U = TRUE, pch.symbols = FALSE)

seq_f1f2_summary <- SummarySeq(SeqList = seq_f1f2)



##### ---- f2_f2 ---- #####

seq_f2 <- sequoia(GenoM = check_thin100K_f2,
                  LifeHistData = LH_f2,
                  Module = "ped",
                  Err = errM,
                  Complex = "full",
                  Herm = "no",
                  UseAge = "yes",
                  args.AP=list(Discrete = TRUE, 
                               MinAgeParent = 2, MaxAgeParent = 2),
                  CalcLLR = TRUE,
                  StrictGenoCheck = TRUE,
                  DummyPrefix = c("F", "M"),
                  Tfilter = -2,
                  Tassign = 1.0)

gmr_f2 <- GetMaybeRel(GenoM = check_thin100K_f2,
                      AgePrior = seq_f2[["AgePriors"]],
                      Err = errM,
                      Module = "ped",
                      Complex = "full",
                      LifeHistData = LH_f2,
                      Herm = "no",
                      quiet = FALSE,
                      Tfilter = -2,
                      Tassign = 1.0,
                      MaxPairs = 7 * nrow(check_thin100K_f2))

relm_f2 <- GetRelM(Pedigree = seq_f2[["Pedigree"]],
                   GenBack = 1, 
                   patmat = FALSE,
                   directed = TRUE,
                   Return = 'Matrix')
table(unique(relm_f2))

relmf2_plot <- PlotRelPairs(RelM = relm_f2, pch.symbols = TRUE)

seq_f2_summary <- SummarySeq(SeqList = seq_f2)

##### ---- f0_f2 ---- #####
seq_f0f2 <- sequoia(GenoM = check_thin100K_f0f2,
                    LifeHistData = LH_f0f2,
                    Module = "ped",
                    Err = errM,
                    Complex = "full",
                    Herm = "no",
                    UseAge = "yes",
                    args.AP=list(Discrete = TRUE, 
                                 MinAgeParent = 2, MaxAgeParent = 2),
                    CalcLLR = TRUE,
                    StrictGenoCheck = TRUE,
                    DummyPrefix = c("F", "M"),
                    Tfilter = -2,
                    Tassign = 1.0)

gmr_f0f2 <- GetMaybeRel(GenoM = check_thin100K_f0f2,
                        AgePrior = seq_f0f2[["AgePriors"]],
                        Err = errM,
                        Module = "ped",
                        Complex = "full",
                        LifeHistData = LH_f0f2,
                        Herm = "no",
                        quiet = FALSE,
                        Tfilter = -2,
                        Tassign = 1.0,
                        MaxPairs = 7 * nrow(check_thin100K_f0f2))

relm_f0f2 <- GetRelM(Pedigree = seq_f0f2[["Pedigree"]],
                     GenBack = 1, 
                     patmat = FALSE,
                     directed = TRUE,
                     Return = 'Matrix')
table(unique(relm_f0f2))

relmf0f2_plot <- PlotRelPairs(RelM = relm_f0f2, drop.U = TRUE, pch.symbols = FALSE)

seq_f0f2_summary <- SummarySeq(SeqList = seq_f0f2)









##### ---- get likelihoods/probs for all pairs of gens ---- #####

# look, we should be using gmr to our advantage. it adds probable assignments that
# didn't get made using sequoia(). i think, given that it's a low div system, and
# that we're UNDERASSIGNING individuals in the test group (assignment rate never
# reached 100%), accuracy was consistently high, so we might have some false positives,
# but still, i think we want to explore all the relationships possible, and use
# probabilities to validate.

# so here's the set of steps.
# 1. make sure that the gmr outputs condition on the pedigrees from above
# then, that output is gonna be ready to compare to the CalcPairLL()/LLtoProb()
# outputs.

# create Pairs df for all six runs
  # start with LH dfs

# run CalcPairLL for all six runs

# run LLtoProb on all six runs
  # all pairwise combs with probs

# validate the toprel for each method for the pairs in the gmr output using a for loop

# if all is clear, plot the probs of all the pairwise ones, plot the assigned ones.

##### ---- cross check the assignments with the TopRel from CalcPairLL() ---- #####

# NOW is the place where you take the assignments and the inferred relationships
# from sequoia() + GetMaybeRel() and you see if that TopRel from the GetMaybeRel()
# output for that pair is the same relationship as the TopRel from CalcPairLL()
# and LLtoProb(). use that as a validation step. plot the pairwise ranges of probs,
# for all inds in that run, and then plot next to it the probs from the inferred
# relationship. 

# i found on 120925 that the LLRs do not match up, between sequoia()/GetMaybeRel()
# (which DO have the same LLRs for a given set of inds), and CalcPairLL, which do
# NOT have the same LLRs for a given set of inds as the other two. that's okay. 
# i can only assume that, if it's not intentional, then it's one of the unavoidable
# reasons why the LLRs are not interpretable. 




























