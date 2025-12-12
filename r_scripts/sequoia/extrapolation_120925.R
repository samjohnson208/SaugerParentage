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
                     Tassign = 0.5)

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
                  Tassign = 0.5)
# assigned 28 dams and 28 sires to 250 + 28 individuals (real + dummy)

gmr_f0 <- GetMaybeRel(GenoM = check_thin100K_f0,
                      SeqList = seq_f0,
                      AgePrior = seq_f0[["AgePriors"]],
                      Err = errM,
                      Module = "ped",
                      Complex = "full",
                      LifeHistData = LH_f0,
                      Herm = "no",
                      quiet = FALSE,
                      Tfilter = -2,
                      Tassign = 0.5,
                      MaxPairs = 7 * nrow(check_thin100K_f0))
# Found 0 likely parent-offspring pairs, and 127, other non-assigned pairs of possible relatives

relm_f0 <- GetRelM(Pedigree = seq_f0[["Pedigree"]],
                   Pairs = gmr_f0$MaybeRel,
                   GenBack = 1, 
                   patmat = FALSE,
                   directed = TRUE,
                   Return = 'Matrix')
table(unique(relm_f0))

relmf0_plot <- PlotRelPairs(RelM = relm_f0, 
                            drop.U = TRUE, 
                            pch.symbols = TRUE,
                            cex.axis = 0.3,
                            mar = c(5, 5, 1, 8))

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
                    Tassign = 0.5)
# assigned 78 dams and 73 sires to 576 + 44 individuals (real + dummy)

gmr_f0f1 <- GetMaybeRel(GenoM = check_thin100K_f0f1,
                        SeqList = seq_f0f1,
                        AgePrior = seq_f0f1[["AgePriors"]],
                        Err = errM,
                        Module = "ped",
                        Complex = "full",
                        LifeHistData = LH_f0f1,
                        Herm = "no",
                        quiet = FALSE,
                        Tfilter = -2,
                        Tassign = 0.5,
                        MaxPairs = 7 * nrow(check_thin100K_f0f1))
# Found 3 likely parent-offspring pairs, and 422, other non-assigned pairs of possible relatives

relm_f0f1 <- GetRelM(Pedigree = seq_f0f1[["Pedigree"]],
                     Pairs = gmr_f0f1$MaybeRel,
                     GenBack = 1, 
                     patmat = FALSE,
                     directed = TRUE,
                     Return = 'Matrix')
table(unique(relm_f0f1))

relmf0f1_plot <- PlotRelPairs(RelM = relm_f0f1, 
                              drop.U = TRUE, 
                              pch.symbols = TRUE,
                              cex.axis = 0.3,
                              mar = c(5, 5, 1, 8))
# no way are these plots useful if there are this many inds on each axis.

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
                  Tassign = 0.5)
# assigned 20 dams and 20 sires to 326 + 18 individuals (real + dummy)

gmr_f1 <- GetMaybeRel(GenoM = check_thin100K_f1,
                      SeqList = seq_f1,
                      AgePrior = seq_f1[["AgePriors"]],
                      Err = errM,
                      Module = "ped",
                      Complex = "full",
                      LifeHistData = LH_f1,
                      Herm = "no",
                      quiet = FALSE,
                      Tfilter = -2,
                      Tassign = 0.5,
                      MaxPairs = 7 * nrow(check_thin100K_f1))
# Found 0 likely parent-offspring pairs, and 170, other non-assigned pairs of possible relatives

relm_f1 <- GetRelM(Pedigree = seq_f1[["Pedigree"]],
                   Pairs = gmr_f1$MaybeRel,
                   GenBack = 1, 
                   patmat = FALSE,
                   directed = TRUE,
                   Return = 'Matrix')
table(unique(relm_f1))

relmf1_plot <- PlotRelPairs(RelM = relm_f1, 
                            drop.U = TRUE,
                            pch.symbols = TRUE,
                            cex.axis = 0.3,
                            mar = c(5, 5, 1, 8))

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
                    Tassign = 0.5)
# assigned 61 dams and 61 sires to 780 + 56 individuals (real + dummy) 

gmr_f1f2 <- GetMaybeRel(GenoM = check_thin100K_f1f2,
                        SeqList = seq_f1f2,
                        AgePrior = seq_f1f2[["AgePriors"]],
                        Err = errM,
                        Module = "ped",
                        Complex = "full",
                        LifeHistData = LH_f1f2,
                        Herm = "no",
                        quiet = FALSE,
                        Tfilter = -2,
                        Tassign = 0.5,
                        MaxPairs = 7 * nrow(check_thin100K_f1f2))
# Found 43 likely parent-offspring pairs, and 552, other non-assigned pairs of possible relatives
# Found 1 parent-parent-offspring trios

relm_f1f2 <- GetRelM(Pedigree = seq_f1f2[["Pedigree"]],
                     Pairs = gmr_f1f2$MaybeRel,
                     GenBack = 1, 
                     patmat = FALSE,
                     directed = TRUE,
                     Return = 'Matrix')
table(unique(relm_f1f2))

relmf1f2_plot <- PlotRelPairs(RelM = relm_f1f2, 
                              drop.U = TRUE, 
                              pch.symbols = TRUE,
                              cex.axis = 0.3,
                              mar = c(5, 5, 1, 8))

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
                               MinAgeParent = 1, MaxAgeParent = 1),
                  CalcLLR = TRUE,
                  StrictGenoCheck = TRUE,
                  DummyPrefix = c("F", "M"),
                  Tfilter = -2,
                  Tassign = 0.5)
# assigned 39 dams and 39 sires to 454 + 36 individuals (real + dummy) 


gmr_f2 <- GetMaybeRel(GenoM = check_thin100K_f2,
                      SeqList = seq_f2,
                      AgePrior = seq_f2[["AgePriors"]],
                      Err = errM,
                      Module = "ped",
                      Complex = "full",
                      LifeHistData = LH_f2,
                      Herm = "no",
                      quiet = FALSE,
                      Tfilter = -2,
                      Tassign = 0.5,
                      MaxPairs = 7 * nrow(check_thin100K_f2))
# Found 0 likely parent-offspring pairs, and 253, other non-assigned pairs of possible relatives


relm_f2 <- GetRelM(Pedigree = seq_f2[["Pedigree"]],
                   Pairs = gmr_f2$MaybeRel,
                   GenBack = 1, 
                   patmat = FALSE,
                   directed = TRUE,
                   Return = 'Matrix')
table(unique(relm_f2))

relmf2_plot <- PlotRelPairs(RelM = relm_f2, 
                            drop.U = TRUE,
                            pch.symbols = TRUE,
                            cex.axis = 0.3,
                            mar = c(5, 5, 1, 8))

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
                    Tassign = 0.5)
# assigned 78 dams and 74 sires to 704 + 65 individuals (real + dummy) 


gmr_f0f2 <- GetMaybeRel(GenoM = check_thin100K_f0f2,
                        SeqList = seq_f0f2,
                        AgePrior = seq_f0f2[["AgePriors"]],
                        Err = errM,
                        Module = "ped",
                        Complex = "full",
                        LifeHistData = LH_f0f2,
                        Herm = "no",
                        quiet = FALSE,
                        Tfilter = -2,
                        Tassign = 0.5,
                        MaxPairs = 7 * nrow(check_thin100K_f0f2))
# Found 1 likely parent-offspring pairs, and 531, other non-assigned pairs of possible relatives

relm_f0f2 <- GetRelM(Pedigree = seq_f0f2[["Pedigree"]],
                     Pairs = gmr_f0f2$MaybeRel,
                     GenBack = 1, 
                     patmat = FALSE,
                     directed = TRUE,
                     Return = 'Matrix')
table(unique(relm_f0f2))

relmf0f2_plot <- PlotRelPairs(RelM = relm_f0f2, 
                              drop.U = TRUE, 
                              pch.symbols = TRUE,
                              cex.axis = 0.3,
                              mar = c(5, 5, 1, 8))

seq_f0f2_summary <- SummarySeq(SeqList = seq_f0f2)
# looks like a few PO duos between f0 and f2. crazy.








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
# outputs. DONE.

      # gmr_f0
      # gmr_f0f1
      # gmr_f1
      # gmr_f1f2
      # gmr_f2
      # gmr_f0f2

# create Pairs df for all six runs
  # start with LH dfs

##### ---- Pairs_f0 ---- #####
# creating the pairs df: create all combinations of id1 and 2, remove rows where 
# they're the same
library(tidyr)
IDs <- rownames(check_thin100K_f0)
length(IDs)
Pairs_f0 <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_f0 <- Pairs_f0 %>% 
  left_join(LH_f0 %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_f0 %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_f0$AgeDif <- Pairs_f0$BY2 - Pairs_f0$BY1

Pairs_f0$focal <- "U"

dim(Pairs_f0)


##### ---- Pairs_f0f1 ---- #####
# creating the pairs df: create all combinations of id1 and 2, remove rows where 
# they're the same
IDs <- rownames(check_thin100K_f0f1)
length(IDs)
Pairs_f0f1 <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_f0f1 <- Pairs_f0f1 %>% 
  left_join(LH_f0f1 %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_f0f1 %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_f0f1$AgeDif <- Pairs_f0f1$BY2 - Pairs_f0f1$BY1

Pairs_f0f1$focal <- "U"

dim(Pairs_f0f1)

##### ---- Pairs_f1 ---- #####
IDs <- rownames(check_thin100K_f1)
length(IDs)
Pairs_f1 <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_f1 <- Pairs_f1 %>% 
  left_join(LH_f1 %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_f1 %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_f1$AgeDif <- Pairs_f1$BY2 - Pairs_f1$BY1

Pairs_f1$focal <- "U"

dim(Pairs_f1)

##### ---- Pairs_f1f2 ---- #####
IDs <- rownames(check_thin100K_f1f2)
length(IDs)
Pairs_f1f2 <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_f1f2 <- Pairs_f1f2 %>% 
  left_join(LH_f1f2 %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_f1f2 %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_f1f2$AgeDif <- Pairs_f1f2$BY2 - Pairs_f1f2$BY1

Pairs_f1f2$focal <- "U"

dim(Pairs_f1f2)

##### ---- Pairs_f2---- #####
IDs <- rownames(check_thin100K_f2)
length(IDs)
Pairs_f2 <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_f2 <- Pairs_f2 %>% 
  left_join(LH_f2 %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_f2 %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_f2$AgeDif <- Pairs_f2$BY2 - Pairs_f2$BY1

Pairs_f2$focal <- "U"

dim(Pairs_f2)

##### ---- Pairs_f0f2 ---- #####
IDs <- rownames(check_thin100K_f0f2)
length(IDs)
Pairs_f0f2 <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_f0f2 <- Pairs_f0f2 %>% 
  left_join(LH_f0f2 %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_f0f2 %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_f0f2$AgeDif <- Pairs_f0f2$BY2 - Pairs_f0f2$BY1

Pairs_f0f2$focal <- "U"

dim(Pairs_f0f2)

##### ---- Pairs_Test ---- #####
IDs <- rownames(check_thin100K_test)
Pairs_test <- expand_grid(
  ID1 = IDs,
  ID2 = IDs) %>%  
  dplyr::filter(ID1 != ID2)

# now join that with the LH_Data so you can get the birth years, sex, 
Pairs_test <- Pairs_test %>% 
  left_join(LH_Test %>% select(ID, Sex, BirthYear),
            by = c("ID1" = "ID")) %>% 
  rename(Sex1 = Sex, BY1 = BirthYear) %>% 
  
  left_join(LH_Test %>% select(ID, Sex, BirthYear),
            by = c("ID2" = "ID")) %>% 
  rename(Sex2 = Sex, BY2 = BirthYear)

Pairs_test$AgeDif <- Pairs_test$BY2 - Pairs_test$BY1

Pairs$focal <- "U"

dim(Pairs_test)
##### ---- ---- #####


##### ---- Getting LLRs and probs for all relationships ---- #####

# run CalcPairLL for all six runs
  # DON'T CONDITION ON A PEDIGREE, THESE ARE ALL POSSIBLE PAIRS
  # THEN WE SAY, HERE ARE THE ONES THAT ACTUALLY GOT PULLED FROM seq, gmr()

##### ---- PairLL_f0 -> prob_pairs_f0  ---- #####

PairLL_f0 <- CalcPairLL(Pairs = Pairs_f0,
                        GenoM = check_thin100K_f0,
                        LifeHistData = LH_f0,
                        AgePrior = seq_f0[["AgePriors"]],
                        Module = "ped",
                        Complex = "full",
                        Herm = 'no',
                        InclDup = FALSE,
                        Err = errM,
                        Tassign = 0.5,
                        Tfilter = -2,
                        quiet = FALSE,
                        Plot = TRUE)
prob_pairs_f0 <- plyr::aaply(as.matrix(PairLL_f0[,10:16]), .margin = 1, LLtoProb)
prob_pairs_f0 <- cbind(PairLL_f0[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_f0)

##### ---- PairLL_f0f1 -> prob_pairs_f0f1 ---- #####

PairLL_f0f1 <- CalcPairLL(Pairs = Pairs_f0f1,
                          GenoM = check_thin100K_f0f1,
                          LifeHistData = LH_f0f1,
                          AgePrior = seq_f0f1[["AgePriors"]],
                          Module = "ped",
                          Complex = "full",
                          Herm = 'no',
                          InclDup = FALSE,
                          Err = errM,
                          Tassign = 0.5,
                          Tfilter = -2,
                          quiet = FALSE,
                          Plot = TRUE)
prob_pairs_f0f1 <- plyr::aaply(as.matrix(PairLL_f0f1[,10:16]), .margin = 1, LLtoProb)
prob_pairs_f0f1 <- cbind(PairLL_f0f1[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_f0f1)

# setwd("/Users/samjohnson/Desktop/")
# save.image(file = "EOD_121025.RData")

##### ---- PairLL_f1 -> prob_pairs_f1 ---- #####

PairLL_f1 <- CalcPairLL(Pairs = Pairs_f1,
                        GenoM = check_thin100K_f1,
                        LifeHistData = LH_f1,
                        AgePrior = seq_f1[["AgePriors"]],
                        Module = "ped",
                        Complex = "full",
                        Herm = 'no',
                        InclDup = FALSE,
                        Err = errM,
                        Tassign = 0.5,
                        Tfilter = -2,
                        quiet = FALSE,
                        Plot = TRUE)
prob_pairs_f1 <- plyr::aaply(as.matrix(PairLL_f1[,10:16]), .margin = 1, LLtoProb)
prob_pairs_f1 <- cbind(PairLL_f1[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_f1)


##### ---- PairLL_f1f2 -> prob_pairs_f1f2 ---- #####

PairLL_f1f2 <- CalcPairLL(Pairs = Pairs_f1f2,
                          GenoM = check_thin100K_f1f2,
                          LifeHistData = LH_f1f2,
                          AgePrior = seq_f1f2[["AgePriors"]],
                          Module = "ped",
                          Complex = "full",
                          Herm = 'no',
                          InclDup = FALSE,
                          Err = errM,
                          Tassign = 0.5,
                          Tfilter = -2,
                          quiet = FALSE,
                          Plot = TRUE)
prob_pairs_f1f2 <- plyr::aaply(as.matrix(PairLL_f1f2[,10:16]), .margin = 1, LLtoProb)
prob_pairs_f1f2 <- cbind(PairLL_f1f2[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_f1f2)

##### ---- PairLL_f2 -> prob_pairs_f2 ---- #####

PairLL_f2 <- CalcPairLL(Pairs = Pairs_f2,
                        GenoM = check_thin100K_f2,
                        LifeHistData = LH_f2,
                        AgePrior = seq_f2[["AgePriors"]],
                        Module = "ped",
                        Complex = "full",
                        Herm = 'no',
                        InclDup = FALSE,
                        Err = errM,
                        Tassign = 0.5,
                        Tfilter = -2,
                        quiet = FALSE,
                        Plot = TRUE)
prob_pairs_f2 <- plyr::aaply(as.matrix(PairLL_f2[,10:16]), .margin = 1, LLtoProb)
prob_pairs_f2 <- cbind(PairLL_f2[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_f2)


# first thing tomorrow. continue with setting up these PairLL_gen dataframes.

# validate the toprel for each method for the pairs in the gmr output using a for loop
  # remember, we're going to have to remove the dummy inds from those relatedness matrices

# answer this question. are the toprels in each the same? are the likelihoods the same?
# pairs represent the rows, have the toprel from each method, have the LLR from each,
# and the prob from method 2, then you want columns that tell you T/F for toprel 
# and LLR being the same between the two methods.

# do the join where you have na's for the pairs that didn't get assigned, and then
# you can subset that dataframe for is.na(TopRel_M1) == FALSE to get the cases
# where the relationships got assigned.

# that way you can work with that whole joined df to generate the only plot you need

# depending on that... decide on what's gonna go on with this stuff below.
# if all is clear, plot the probs of all the pairwise ones, plot the assigned ones.
  # remember, going to have to filter out the agedif = 0 pairs from the prob_pairs
  # dfs that span generations, since we'll have already plotted those. i think?



##### ---- PairLL_f0f2 -> prob_pairs_f0f2 ---- #####

PairLL_f0f2 <- CalcPairLL(Pairs = Pairs_f0f2,
                          GenoM = check_thin100K_f0f2,
                          LifeHistData = LH_f0f2,
                          AgePrior = seq_f0f2[["AgePriors"]],
                          Module = "ped",
                          Complex = "full",
                          Herm = 'no',
                          InclDup = FALSE,
                          Err = errM,
                          Tassign = 0.5,
                          Tfilter = -2,
                          quiet = FALSE,
                          Plot = TRUE)
prob_pairs_f0f2 <- plyr::aaply(as.matrix(PairLL_f0f2[,10:16]), .margin = 1, LLtoProb)
prob_pairs_f0f2 <- cbind(PairLL_f0f2[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_f0f2)

##### ---- PairLL_test <- prob_pairs_test ---- #####


# IMPORTANT NOTE: I did the test group AFTER all the rest of the combinations in
# this script, seq_test, gmr_test, and relm_test were all made below in the "a quick 
# check on the thresholds and test group" section

PairLL_test <- CalcPairLL(Pairs = Pairs_test,
                          GenoM = check_thin100K_test,
                          LifeHistData = LH_Test,
                          AgePrior = seq_test[["AgePriors"]],
                          Module = "ped",
                          Complex = "full",
                          Herm = 'no',
                          InclDup = FALSE,
                          Err = errM,
                          Tassign = 0.5,
                          Tfilter = -2,
                          quiet = FALSE,
                          Plot = TRUE)
prob_pairs_test <- plyr::aaply(as.matrix(PairLL_test[,10:16]), .margin = 1, LLtoProb)
prob_pairs_test <- cbind(PairLL_test[, c("ID1", "ID2","AgeDif", "TopRel")], prob_pairs_test)


##### ---- after creating those... ---- #####

setwd("/Users/samjohnson/Desktop/")
save.image(file = "M1M2_completed_121125.RData")


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

# so we need to create a df that has the following columns

# ind1 ind2 pair group_ind1 group_ind2 group_pair TopRel_seq TopRel_pairwise TopRel_agree Prob_pairwise

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

##### ---- generating results: f0f1 ---- ###### 

# step 1: remove dummy individuals from the relatedness matrix

true_inds <- rownames(relm_f0f1)[str_starts(rownames(relm_f0f1), "SAR")]
length(true_inds) # make sure it's the correct length
relm_f0f1_nodum <- relm_f0f1[true_inds, true_inds]

# step 1.5: we need to clean the relatedness matrix so that when we merge, we can
# get the probabilities for the relationships predicted by method 1
clean_seq_matrix <- function(mat) {
  mat <- as.matrix(mat)  # ensure it's a matrix
  mat[mat %in% c("FA?", "FA")] <- "FA"
  mat[mat %in% c("FS?", "FS")] <- "FS"
  mat[mat %in% c("HS?", "HS")] <- "HS"
  mat[mat %in% c("O", "PO?", "PO")] <- "PO"
  mat[mat %in% c("Q?", "Q")] <- "??"
  return(mat)
}

relm_f0f1_clean <- clean_seq_matrix(relm_f0f1_nodum)

# step 2: convert the relatedness matrix to a dataframe, long format with unique, sorted pairs

relm_f0f1_long <- relm_f0f1_clean %>%
    as.data.frame() %>% # turn to df
    tibble::rownames_to_column("ID1") %>% # set up the first col as id1
    pivot_longer(cols = -ID1, # then pivot so you set up everything else around id1
                 names_to = "ID2",
                 values_to = "TopRel_seq") %>% 
  # now you have the pairs and their relationship from this method in 3 columns
  # remove self comparisons
   filter(ID1 != ID2) %>%
  # create sorted pair IDs so duplicates collapse
   mutate(ind1 = pmin(ID1, ID2),
          ind2 = pmax(ID1, ID2)) %>%
   distinct(ind1, ind2, .keep_all = TRUE) %>% # make sure you keep the relationship col
   select(ind1, ind2, TopRel_seq) # and place them in this order, if they're not already

# step 3: prepare the prob_pairs df for merging. keep only distinct, sorted pairs

prob_long <- prob_pairs_f0f1 %>%
  # create sorted pair IDs so it matches the matrix, same scheme as above
  mutate(ind1 = pmin(ID1, ID2),
         ind2 = pmax(ID1, ID2)) %>%
  distinct(ind1, ind2, .keep_all = TRUE) # make sure to keep everything else

  # check your dims to make sure these are ready for the merge
  dim(prob_long) 
  dim(relm_f0f1_long)

# step 4: perform the merge

merged <- relm_f0f1_long %>%
  left_join(prob_long,
            by = c("ind1", "ind2")) # join them by each unique pair

# step 5: add the generation to each individual

# set up the assign_gen function for use later
assign_gen <- function(x) {
  case_when(
    x %in% f0_inds$ID ~ "F0",
    x %in% f1_inds$ID ~ "F1",
    x %in% f2_inds$ID ~ "F2",
    TRUE ~ NA_character_ # shouldn't be the case, but put an NA if an ind isn't there
  )
}

# now run assign_gen on all inds in the merged df
merged2 <- merged %>%
  mutate(group_ind1 = assign_gen(ind1),
         group_ind2 = assign_gen(ind2),
         # now we'll use our same pmin and pmax setup to create a group_pair col
         group_pair = paste(pmin(group_ind1, group_ind2),
                            pmax(group_ind1, group_ind2),
                            sep = "_"))

# step 6: do the TopRel's agree? extract the probability of the TopRel_pairwise,
        # and extract the probability of the TopRel_seq
merged3 <- merged2 %>%
  mutate(TopRel_pairwise = TopRel,      # rename the pairwise one for clarity first
         TopRel_agree = TopRel_seq == TopRel_pairwise) %>% # do the TopRel's agree? 
  # now go row by row, pull the TopRel from each method, and give me it's probability
  # in columns that correspond to each method
  rowwise() %>%
  mutate(Prob_pairwise = case_when(
         TopRel_pairwise == "PO" ~ PO,
         TopRel_pairwise == "FS" ~ FS,
         TopRel_pairwise == "HS" ~ HS,
         TopRel_pairwise == "GP" ~ GP,
         TopRel_pairwise == "FA" ~ FA,
         TopRel_pairwise == "HA" ~ HA,
         TopRel_pairwise == "U"  ~ U,
         TRUE ~ NA_real_)) %>%
  mutate(Prob_seq = case_when(
         TopRel_seq == "PO" ~ PO,
         TopRel_seq == "FS" ~ FS,
         TopRel_seq == "HS" ~ HS,
         TopRel_seq == "GP" ~ GP,
         TopRel_seq == "FA" ~ FA,
         TopRel_seq == "HA" ~ HA,
         TopRel_seq == "U"  ~ U,
         TRUE ~ NA_real_)) %>% 
  ungroup()

# step 7: reorder the columns to make it visually navigable

resultstoplot_f0f1 <- merged3 %>%
    select(ind1, ind2, group_ind1, group_ind2, group_pair, TopRel_agree, 
           TopRel_seq, TopRel_pairwise, Prob_seq, Prob_pairwise, everything()) %>% 
    select(-ID1, -ID2, -TopRel)

dim(resultstoplot_f0f1) # double check the dim

##### ---- ---- #####

##### ---- generating results: f1f2 ---- #####

# step 1: remove dummy individuals from the relatedness matrix

true_inds <- rownames(relm_f1f2)[str_starts(rownames(relm_f1f2), "SAR")]
length(true_inds) # make sure it's the correct length
relm_f1f2_nodum <- relm_f1f2[true_inds, true_inds]

# step 1.5: we need to clean the relatedness matrix so that when we merge, we can
# get the probabilities for the relationships predicted by method 1

relm_f1f2_clean <- clean_seq_matrix(relm_f1f2_nodum)

# step 2: convert the relatedness matrix to a dataframe, long format with unique, sorted pairs

relm_f1f2_long <- relm_f1f2_clean %>%
  as.data.frame() %>% # turn to df
  tibble::rownames_to_column("ID1") %>% # set up the first col as id1
  pivot_longer(cols = -ID1, # then pivot so you set up everything else around id1
               names_to = "ID2",
               values_to = "TopRel_seq") %>% 
  # now you have the pairs and their relationship from this method in 3 columns
  # remove self comparisons
  filter(ID1 != ID2) %>%
  # create sorted pair IDs so duplicates collapse
  mutate(ind1 = pmin(ID1, ID2),
         ind2 = pmax(ID1, ID2)) %>%
  distinct(ind1, ind2, .keep_all = TRUE) %>% # make sure you keep the relationship col
  select(ind1, ind2, TopRel_seq) # and place them in this order, if they're not already

# step 3: prepare the prob_pairs df for merging. keep only distinct, sorted pairs

prob_long <- prob_pairs_f1f2 %>%
  # create sorted pair IDs so it matches the matrix, same scheme as above
  mutate(ind1 = pmin(ID1, ID2),
         ind2 = pmax(ID1, ID2)) %>%
  distinct(ind1, ind2, .keep_all = TRUE) # make sure to keep everything else

# check your dims to make sure these are ready for the merge
dim(prob_long) 
dim(relm_f1f2_long)

# step 4: perform the merge

merged <- relm_f1f2_long %>%
  left_join(prob_long,
            by = c("ind1", "ind2")) # join them by each unique pair

# step 5: add the generation to each individual

# set up the assign_gen function for use later
assign_gen <- function(x) {
  case_when(
    x %in% f0_inds$ID ~ "F0",
    x %in% f1_inds$ID ~ "F1",
    x %in% f2_inds$ID ~ "F2",
    TRUE ~ NA_character_ # shouldn't be the case, but put an NA if an ind isn't there
  )
}

# now run assign_gen on all inds in the merged df
merged2 <- merged %>%
  mutate(group_ind1 = assign_gen(ind1),
         group_ind2 = assign_gen(ind2),
         # now we'll use our same pmin and pmax setup to create a group_pair col
         group_pair = paste(pmin(group_ind1, group_ind2),
                            pmax(group_ind1, group_ind2),
                            sep = "_"))

# step 6: do the TopRel's agree? extract the probability of the TopRel_pairwise,
# and extract the probability of the TopRel_seq
merged3 <- merged2 %>%
  mutate(TopRel_pairwise = TopRel,      # rename the pairwise one for clarity first
         TopRel_agree = TopRel_seq == TopRel_pairwise) %>% # do the TopRel's agree? 
  # now go row by row, pull the TopRel from each method, and give me it's probability
  # in columns that correspond to each method
  rowwise() %>%
  mutate(Prob_pairwise = case_when(
    TopRel_pairwise == "PO" ~ PO,
    TopRel_pairwise == "FS" ~ FS,
    TopRel_pairwise == "HS" ~ HS,
    TopRel_pairwise == "GP" ~ GP,
    TopRel_pairwise == "FA" ~ FA,
    TopRel_pairwise == "HA" ~ HA,
    TopRel_pairwise == "U"  ~ U,
    TRUE ~ NA_real_)) %>%
  mutate(Prob_seq = case_when(
    TopRel_seq == "PO" ~ PO,
    TopRel_seq == "FS" ~ FS,
    TopRel_seq == "HS" ~ HS,
    TopRel_seq == "GP" ~ GP,
    TopRel_seq == "FA" ~ FA,
    TopRel_seq == "HA" ~ HA,
    TopRel_seq == "U"  ~ U,
    TRUE ~ NA_real_)) %>% 
  ungroup()

# step 7: reorder the columns to make it visually navigable

resultstoplot_f1f2 <- merged3 %>%
  select(ind1, ind2, group_ind1, group_ind2, group_pair, TopRel_agree, 
         TopRel_seq, TopRel_pairwise, Prob_seq, Prob_pairwise, everything()) %>% 
  select(-ID1, -ID2, -TopRel)

dim(resultstoplot_f1f2) # double check the dim
##### ---- ---- #####

##### ---- generating results: f0f2 ---- #####

# step 1: remove dummy individuals from the relatedness matrix

true_inds <- rownames(relm_f0f2)[str_starts(rownames(relm_f0f2), "SAR")]
length(true_inds) # make sure it's the correct length
relm_f0f2_nodum <- relm_f0f2[true_inds, true_inds]

# step 1.5: we need to clean the relatedness matrix so that when we merge, we can
# get the probabilities for the relationships predicted by method 1

relm_f0f2_clean <- clean_seq_matrix(relm_f0f2_nodum)

# step 2: convert the relatedness matrix to a dataframe, long format with unique, sorted pairs

relm_f0f2_long <- relm_f0f2_clean %>%
  as.data.frame() %>% # turn to df
  tibble::rownames_to_column("ID1") %>% # set up the first col as id1
  pivot_longer(cols = -ID1, # then pivot so you set up everything else around id1
               names_to = "ID2",
               values_to = "TopRel_seq") %>% 
  # now you have the pairs and their relationship from this method in 3 columns
  # remove self comparisons
  filter(ID1 != ID2) %>%
  # create sorted pair IDs so duplicates collapse
  mutate(ind1 = pmin(ID1, ID2),
         ind2 = pmax(ID1, ID2)) %>%
  distinct(ind1, ind2, .keep_all = TRUE) %>% # make sure you keep the relationship col
  select(ind1, ind2, TopRel_seq) # and place them in this order, if they're not already

# step 3: prepare the prob_pairs df for merging. keep only distinct, sorted pairs

prob_long <- prob_pairs_f0f2 %>%
  # create sorted pair IDs so it matches the matrix, same scheme as above
  mutate(ind1 = pmin(ID1, ID2),
         ind2 = pmax(ID1, ID2)) %>%
  distinct(ind1, ind2, .keep_all = TRUE) # make sure to keep everything else

# check your dims to make sure these are ready for the merge
dim(prob_long) 
dim(relm_f0f2_long)

# step 4: perform the merge

merged <- relm_f0f2_long %>%
  left_join(prob_long,
            by = c("ind1", "ind2")) # join them by each unique pair

# step 5: add the generation to each individual

# set up the assign_gen function for use later
assign_gen <- function(x) {
  case_when(
    x %in% f0_inds$ID ~ "F0",
    x %in% f1_inds$ID ~ "F1",
    x %in% f2_inds$ID ~ "F2",
    TRUE ~ NA_character_ # shouldn't be the case, but put an NA if an ind isn't there
  )
}

# now run assign_gen on all inds in the merged df
merged2 <- merged %>%
  mutate(group_ind1 = assign_gen(ind1),
         group_ind2 = assign_gen(ind2),
         # now we'll use our same pmin and pmax setup to create a group_pair col
         group_pair = paste(pmin(group_ind1, group_ind2),
                            pmax(group_ind1, group_ind2),
                            sep = "_"))

# step 6: do the TopRel's agree? extract the probability of the TopRel_pairwise,
# and extract the probability of the TopRel_seq
merged3 <- merged2 %>%
  mutate(TopRel_pairwise = TopRel,      # rename the pairwise one for clarity first
         TopRel_agree = TopRel_seq == TopRel_pairwise) %>% # do the TopRel's agree? 
  # now go row by row, pull the TopRel from each method, and give me it's probability
  # in columns that correspond to each method
  rowwise() %>%
  mutate(Prob_pairwise = case_when(
    TopRel_pairwise == "PO" ~ PO,
    TopRel_pairwise == "FS" ~ FS,
    TopRel_pairwise == "HS" ~ HS,
    TopRel_pairwise == "GP" ~ GP,
    TopRel_pairwise == "FA" ~ FA,
    TopRel_pairwise == "HA" ~ HA,
    TopRel_pairwise == "U"  ~ U,
    TRUE ~ NA_real_)) %>%
  mutate(Prob_seq = case_when(
    TopRel_seq == "PO" ~ PO,
    TopRel_seq == "FS" ~ FS,
    TopRel_seq == "HS" ~ HS,
    TopRel_seq == "GP" ~ GP,
    TopRel_seq == "FA" ~ FA,
    TopRel_seq == "HA" ~ HA,
    TopRel_seq == "U"  ~ U,
    TRUE ~ NA_real_)) %>% 
  ungroup()

# step 7: reorder the columns to make it visually navigable

resultstoplot_f0f2 <- merged3 %>%
  select(ind1, ind2, group_ind1, group_ind2, group_pair, TopRel_agree, 
         TopRel_seq, TopRel_pairwise, Prob_seq, Prob_pairwise, everything()) %>% 
  select(-ID1, -ID2, -TopRel)

dim(resultstoplot_f0f2) # double check the dim


##### ---- ---- #####

##### ---- generating results: test ---- #####

# step 1: remove dummy individuals from the relatedness matrix

true_inds <- rownames(relm_test)[str_starts(rownames(relm_test), "SAR")]
length(true_inds) # make sure it's the correct length
relm_test_nodum <- relm_test[true_inds, true_inds]

# step 1.5: we need to clean the relatedness matrix so that when we merge, we can
# get the probabilities for the relationships predicted by method 1

relm_test_clean <- clean_seq_matrix(relm_test_nodum)

# step 2: convert the relatedness matrix to a dataframe, long format with unique, sorted pairs

relm_test_long <- relm_test_clean %>%
  as.data.frame() %>% # turn to df
  tibble::rownames_to_column("ID1") %>% # set up the first col as id1
  pivot_longer(cols = -ID1, # then pivot so you set up everything else around id1
               names_to = "ID2",
               values_to = "TopRel_seq") %>% 
  # now you have the pairs and their relationship from this method in 3 columns
  # remove self comparisons
  filter(ID1 != ID2) %>%
  # create sorted pair IDs so duplicates collapse
  mutate(ind1 = pmin(ID1, ID2),
         ind2 = pmax(ID1, ID2)) %>%
  distinct(ind1, ind2, .keep_all = TRUE) %>% # make sure you keep the relationship col
  select(ind1, ind2, TopRel_seq) # and place them in this order, if they're not already

# step 3: prepare the prob_pairs df for merging. keep only distinct, sorted pairs

prob_long <- prob_pairs_test %>%
  # create sorted pair IDs so it matches the matrix, same scheme as above
  mutate(ind1 = pmin(ID1, ID2),
         ind2 = pmax(ID1, ID2)) %>%
  distinct(ind1, ind2, .keep_all = TRUE) # make sure to keep everything else

# check your dims to make sure these are ready for the merge
dim(prob_long) 
dim(relm_test_long)

# step 4: perform the merge

merged <- relm_test_long %>%
  left_join(prob_long,
            by = c("ind1", "ind2")) # join them by each unique pair

# step 5: add the generation to each individual

# set up the assign_gen function for use later
assign_gen <- function(x) {
  case_when(
    x %in% f0_inds$ID ~ "F0",
    x %in% f1_inds$ID ~ "F1",
    x %in% f2_inds$ID ~ "F2",
    TRUE ~ NA_character_ # shouldn't be the case, but put an NA if an ind isn't there
  )
}

# now run assign_gen on all inds in the merged df
merged2 <- merged %>%
  mutate(group_ind1 = assign_gen(ind1),
         group_ind2 = assign_gen(ind2),
         # now we'll use our same pmin and pmax setup to create a group_pair col
         group_pair = paste(pmin(group_ind1, group_ind2),
                            pmax(group_ind1, group_ind2),
                            sep = "_"))

# step 6: do the TopRel's agree? extract the probability of the TopRel_pairwise,
# and extract the probability of the TopRel_seq
merged3 <- merged2 %>%
  mutate(TopRel_pairwise = TopRel,      # rename the pairwise one for clarity first
         TopRel_agree = TopRel_seq == TopRel_pairwise) %>% # do the TopRel's agree? 
  # now go row by row, pull the TopRel from each method, and give me it's probability
  # in columns that correspond to each method
  rowwise() %>%
  mutate(Prob_pairwise = case_when(
    TopRel_pairwise == "PO" ~ PO,
    TopRel_pairwise == "FS" ~ FS,
    TopRel_pairwise == "HS" ~ HS,
    TopRel_pairwise == "GP" ~ GP,
    TopRel_pairwise == "FA" ~ FA,
    TopRel_pairwise == "HA" ~ HA,
    TopRel_pairwise == "U"  ~ U,
    TRUE ~ NA_real_)) %>%
  mutate(Prob_seq = case_when(
    TopRel_seq == "PO" ~ PO,
    TopRel_seq == "FS" ~ FS,
    TopRel_seq == "HS" ~ HS,
    TopRel_seq == "GP" ~ GP,
    TopRel_seq == "FA" ~ FA,
    TopRel_seq == "HA" ~ HA,
    TopRel_seq == "U"  ~ U,
    TRUE ~ NA_real_)) %>% 
  ungroup()

# step 7: reorder the columns to make it visually navigable

resultstoplot_test <- merged3 %>%
  select(ind1, ind2, group_ind1, group_ind2, group_pair, TopRel_agree, 
         TopRel_seq, TopRel_pairwise, Prob_seq, Prob_pairwise, everything()) %>% 
  select(-ID1, -ID2, -TopRel)

dim(resultstoplot_test) # double check the dim

##### ---- ---- #####

setwd("/Users/samjohnson/Desktop/")
save.image(file = "resultstoplot_121125.RData")


##### ---- exploring the output ---- #####
# how many times do the relationships disagree/agree among the methods?
table(resultstoplot_f0f1$TopRel_agree)
table(resultstoplot_f1f2$TopRel_agree)
table(resultstoplot_f0f2$TopRel_agree)
table(resultstoplot_test$TopRel_agree)

# how many times do disagreements happen because the GetMaybeRel() relationships have "?" at the end?
table(resultstoplot_f0f1$TopRel_seq)
table(resultstoplot_f0f1$TopRel_pairwise)

table(resultstoplot_f1f2$TopRel_seq)
table(resultstoplot_f1f2$TopRel_pairwise)

table(resultstoplot_f0f2$TopRel_seq)
table(resultstoplot_f0f2$TopRel_pairwise)

table(resultstoplot_test$TopRel_seq)
table(resultstoplot_test$TopRel_pairwise)
##### ---- ---- #####

##### ---- removing the U/U scenario ---- #####

# there are TONS of cases in each dataframe where both methods agree that the inds
# are unrelated, and those clusters of points are clogging up my plots.

resultstoplot_f0f1_cleanUU <- resultstoplot_f0f1 %>% 
    filter(!(TopRel_seq == "U" & TopRel_pairwise == "U"))

resultstoplot_f1f2_cleanUU <- resultstoplot_f1f2 %>% 
  filter(!(TopRel_seq == "U" & TopRel_pairwise == "U"))

resultstoplot_f0f2_cleanUU <- resultstoplot_f0f2 %>% 
  filter(!(TopRel_seq == "U" & TopRel_pairwise == "U"))

table(resultstoplot_f0f1$TopRel_seq)
table(resultstoplot_f0f1$TopRel_pairwise)

##### ---- ---- #####
# ^ I AM TEMPORARILY BYPASSING THIS TO SEE IF I CAN MAKE THE BOXPLOTS WITH THE U/U
# SCENARIOS. I THINK REMOVING THEM WOULD BE MISLEADING.


##### ---- removing the question marks after inferred relationships ---- #####

# sequoia() introduces a question mark after an inferred relationship if it's low
# confidence, and that's making the plots a bit weird. since there are so many points,
# i can't seem to get those situations to plot, and i'd like the relationships and
# probabilities to be transferable between methods. if the methods don't agree, they
# are likely low probability, and i'd like those situations to contribute to the 
# distributions that i'm plotting.

# for example, there are some situations where TopRel_seq == "PO?" and TopRel_pairwise == "PO"
# and in those cases, the Prob_seq is NA, since CalcPairLL and LLtoProb don't make
# probability columns for these relationships with question marks on the end.

# this function and the following code will take care of those situations

clean_seq_relationships <- function(df) {
  df %>%
    mutate(
      TopRel_seq = case_when(
        TopRel_seq %in% c("FA?", "FA") ~ "FA",
        TopRel_seq %in% c("FS?", "FS") ~ "FS",
        TopRel_seq %in% c("HS?", "HS") ~ "HS",
        TopRel_seq %in% c("O", "PO?","PO") ~ "PO",
        TopRel_seq %in% c("Q?") ~ "??",
        TRUE ~ TopRel_seq
      )
    )
}

toplot_f0f1_cleanUU_noQM <- clean_seq_relationships(resultstoplot_f0f1_cleanUU)
toplot_f1f2_cleanUU_noQM <- clean_seq_relationships(resultstoplot_f1f2_cleanUU)
toplot_f0f2_cleanUU_noQM <- clean_seq_relationships(resultstoplot_f0f2_cleanUU)

# didn't realize the "??" originated from the pairwise method. now just need to
# filter for cases where NEITHER of the two methods generate a pair with "??"

##### ---- ---- #####
# ^ I AM ALSO BYPASSING THIS, SINCE I MOVED THIS FUNCTIONALITY TO STEP 1.5 ABOVE


##### ---- remove the unknowns ("??")---- ####
remove_unknowns <- function(df) {
  df %>%
    filter(
      TopRel_seq != "??",
      TopRel_pairwise != "??"
    )
}

toplot_f0f1_noQM <- remove_unknowns(resultstoplot_f0f1)
table(toplot_f0f1_noQM$TopRel_pairwise)
table(toplot_f0f1_noQM$TopRel_seq)

toplot_f1f2_noQM <- remove_unknowns(resultstoplot_f1f2)
table(toplot_f1f2_noQM$TopRel_pairwise)
table(toplot_f1f2_noQM$TopRel_seq)

toplot_f0f2_noQM <- remove_unknowns(resultstoplot_f0f2)
table(toplot_f0f2_noQM$TopRel_pairwise)
table(toplot_f0f2_noQM$TopRel_seq)

toplot_test_noQM <- remove_unknowns(resultstoplot_test)
table(toplot_test_noQM$TopRel_pairwise)
table(toplot_test_noQM$TopRel_seq)

##### ---- ---- #####
# ^ HOWEVER, WE DO STILL NEED TO DO THIS TO GET RID OF THE "??" SCENARIO

##### ---- pivot and prep the resultstoplot_ ---- #####
# we've got to pivot longer for the facet wrap in ggplot to work, so here's a repeatable
# function that we can use for all three of the dataframes

prep_plot_data <- function(df) {
  df %>%
    select(ind1, ind2, TopRel_pairwise, TopRel_seq, Prob_pairwise, Prob_seq) %>% 
    # keep only these few columns
    tidyr::pivot_longer(cols = c(Prob_pairwise, Prob_seq),
                        names_to = "Method",
                        values_to = "Probability") %>% 
    # ^ now we take the values in the probability columns and place them in
    # a single probability column, and you take the names of the columns that we're
    # combining, and we put those in the new "Method" column
    mutate(Method = dplyr::recode(Method,
                             Prob_pairwise = "Pairwise",
                             Prob_seq      = "Sequoia"),
      # ^ now we take that method column, and use this recode function which renames
      # the values in that Method column using this kinda "casewhen" system. Makes
      # the method names readable for the facet wrap.
            Relationship = ifelse(Method == "Pairwise",
                                  TopRel_pairwise,
                                  TopRel_seq)
    # ^ lastly, we need to place the relationships in a corresponding "Relationship"
    # column, so if the method's pairwise, put the TopRel_pairwise in there. if
    # it ain't, put the sequoia in there, which will correspond to Method == "Sequoia"
    )
}

# pivot the resultstoplot_ dataframes, they should double in length
toplot_f0f1_noQM_piv <- prep_plot_data(toplot_f0f1_noQM)
toplot_f1f2_noQM_piv <- prep_plot_data(toplot_f1f2_noQM)
toplot_f0f2_noQM_piv <- prep_plot_data(toplot_f0f2_noQM)
toplot_test_noQM_piv <- prep_plot_data(toplot_test_noQM)


# okay so it's in the correct format now, but the problem is that EARLIER in the
# big 7 step process there were still PO? and stuff in the relatedness matrices
# and ?? in the pairwise dataframes. need to go do these filtering steps earlier
# in the pipeline.
# ^ THIS HAS BEEN RECONCILED, ANYTHING IN THE "plotting results" SECTION IS GOOD
##### ---- ---- #####

# pre-plot data save
setwd("/Users/samjohnson/Desktop/")
save.image(file = "resultstoplot_121125.RData")

##### ---- plotting results ---- ##### 

plot_rel_boxpoints <- function(plotdata, title = NULL) {

  # define the levels of the relationships to order them on the x-axis
  rel_levels <- c("PO", "FA", "HA", "FS", "HS", "??", "U")
  plotdata <- plotdata %>%
    mutate(Relationship = factor(Relationship, levels = rel_levels))
  
  # plot the boxplot and jittered points
  ggplot(plotdata, aes(x = Relationship, y = Probability, fill = Relationship)) +
    
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Hide outliers since we add points
    
    geom_jitter(aes(color = Relationship),
                width = 0.15, height = 0,
                size = 0.9, alpha = 0.5) +
    
    facet_wrap(~ Method, nrow = 1) +
    
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    
    theme_bw() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    
    labs(
      x = "Relationship Category",
      y = "Probability",
      title = title
    )
}

plot_f0f1 <- plot_rel_boxpoints(toplot_f0f1_noQM_piv, "F0–F1: Relationship Probabilities by Method")
plot_f0f1

plot_f1f2 <- plot_rel_boxpoints(toplot_f1f2_noQM_piv, "F1–F2: Relationship Probabilities by Method")
plot_f1f2

plot_f0f2 <- plot_rel_boxpoints(toplot_f0f2_noQM_piv, "F0–F2: Relationship Probabilities by Method")
plot_f0f2

plot_test <- plot_rel_boxpoints(toplot_test_noQM_piv, "Test Group: Relationship Probabilities by Method")
plot_test
##### ---- ---- #####

# post-plotting data save
setwd("/Users/samjohnson/Desktop/")
save.image(file = "resultsplotted_EOD_121125.RData")

##### ----  unsuccessful plotting functions and plot types ---- #####

plot_rel_violin <- function(plotdata, title = NULL) {
  
  # Ensure Relationship is a factor for consistent coloring
  plotdata <- plotdata %>%
    mutate(Relationship = factor(Relationship,
                                 levels = c("??", "FA", "FA?", "FS", "FS?", "HA", "HS", "HS?", "O", "PO", "PO?", "Q?", "U")))
  
  ggplot(plotdata, aes(x = Relationship, y = Probability, fill = Relationship)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = "black") +  # outline for visibility
    # geom_jitter(aes(color = Relationship),
    #             width = 0.15, height = 0,
    #             alpha = 0.8, size = 0.8) +
    facet_wrap(~ Method, nrow = 1) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    theme_bw() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      x = "Top Relationship",
      y = "Probability",
      title = title
    )
}

plot_rel_violin_large <- function(plotdata, title = NULL, sample_frac = 0.02) {
  
  plotdata <- plotdata %>%
    mutate(Relationship = factor(Relationship,
                                 levels = c("??", "FA", "FA?", "FS", "FS?", "HA", "HS", "HS?", "O", "PO", "PO?", "Q?", "U")))
  
  ggplot(plotdata, aes(x = Relationship, y = Probability, fill = Relationship)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = "black") +
    geom_boxplot(width = 0.15, outlier.size = 0.5) +
    geom_jitter(
      data = plotdata %>% sample_frac(sample_frac),  # plot only a fraction of points
      aes(color = Relationship),
      width = 0.15, height = 0,
      alpha = 0.7, size = 0.6
    ) +
    facet_wrap(~ Method, nrow = 1) +
    scale_fill_viridis_d(option = "viridis") +
    scale_color_viridis_d(option = "viridis") +
    theme_bw() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      x = "Top Relationship",
      y = "Probability",
      title = title
    )
}

plotf0f1 <- plot_rel_violin_large(toplot_f0f1_cleanUU_prepped, title = "F0-F1 Relationships and Probabilities by Method")
plotf0f1
dim(toplot_f0f1_cleanUU_prepped)

plotf1f2 <- plot_rel_violin_large(toplot_f1f2_cleanUU_prepped, title = "F1-F2 Relationships and Probabilities by Method")
plotf1f2
dim(toplot_f1f2_cleanUU_prepped)

plotf0f2 <- plot_rel_violin_large(toplot_f0f2_cleanUU_prepped, title = "F0-F2 Relationships and Probabilities by Method")
plotf0f2
dim(toplot_f0f2_cleanUU_prepped)
##### ---- ---- #####

##### ---- a quick check on the thresholds and test group ---- ##### 

# alright, as a refresher, the Tfilter is the first step, which allows relationships
# to even be considered or not. a STRICTER Tfilter is MORE negative than the default,
# which is -2. Tassign decides BETWEEN relationships for a given individual (e.g., 
# between two closely related candidate parents). this value needs to be positive.
# a STRICTER Tassign is MORE positive, (i.e., there needs to be a larger difference 
# in the likelihoods between two relationsihps for one of them to be chosen). that
# is to say, a lower Tassign is LESS conservative/strict.

# why any of this is the case? i don't know. i don't know how these are calculated,
# i don't know how the "ratio" comes into play, since it seems like these thresholds
# are differences rather than ratios. I'm not sure how this filtering process works
# when deciding relationships that get displayed by either sequoia() OR GetMaybeRel()

seq_test <- sequoia(GenoM = check_thin100K_test,
                    LifeHistData = LH_Test,
                    args.AP=list(Discrete = TRUE),
                    UseAge = "yes",
                    Module = "ped",
                    Complex = "full",
                    Err = errM,
                    Herm = "no",
                    CalcLLR = TRUE,
                    StrictGenoCheck = TRUE,
                    DummyPrefix = c("F", "M"),
                    Tfilter = -2,
                    Tassign = 0.5)
# Tfilter = -2, -3, -4, -5, Tassign = 0.5
# assigned 98 dams and 98 sires to 206 + 6 individuals (real + dummy)

# Tfilter = -2, -3, -4, -5, Tassign = 1.0
# assigned 94 dams and 95 sires to 206 + 3 individuals (real + dummy)

# so it just adds a few more dummy inds. nothing drastic.

gmr_test <- GetMaybeRel(GenoM = check_thin100K_test,                    
                        LifeHistData = LH_Test,
                        SeqList = seq_test,
                        Err = errM,                                   
                        Module = "ped",
                        Complex = "full",
                        AgePrior = seq_test[["AgePriors"]],
                        quiet = FALSE,
                        Tassign = 0.5,
                        Tfilter = -2,
                        Herm = "no",
                        MaxPairs = 7*nrow(check_thin100K_test)) 

# assignment and filtering threshold changes don't do anything to the output.
# also ran gmr_test with Tfilter as -2, Tassign as 0.5 and 1.0. does not change
# the amount of relationships that are output.

relm_test <- GetRelM(Pedigree = seq_test[["Pedigree"]],
                     Pairs = gmr_test$MaybeRel,
                     GenBack = 1, 
                     patmat = FALSE,
                     directed = TRUE,
                     Return = 'Matrix')
table(unique(relm_test))

relmtest_plot <- PlotRelPairs(RelM = relm_test,
                              pch.symbols = TRUE,
                              mar = c(5, 5, 1, 8))

seq_test_summary <- SummarySeq(SeqList = seq_test)




##### ---- ---- #####






















