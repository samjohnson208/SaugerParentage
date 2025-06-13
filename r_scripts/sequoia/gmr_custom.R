## sequoia_custom.R by SPJ 061225
## PURPOSE: to load a custom version of sequoia::GetMaybeRel() to alter assignment/mismatch thresholds
## USAGE: Rscript sequoia_custom.R

GetMaybeRel_Custom <- function (GenoM = NULL, SeqList = NULL, Pedigree = NULL, LifeHistData = NULL, 
                                AgePrior = NULL, Module = "par", Complex = "full", Herm = "no", 
                                Err = 1e-04, ErrFlavour = "version2.0", Tassign = 0.01, Tfilter = -10, 
                                MaxPairs = 7 * nrow(GenoM), quiet = FALSE, ParSib = NULL, 
                                MaxMismatch = NA) 
{
  on.exit(.Fortran(sequoia:::deallocall), add = TRUE)
  if (!is.null(ParSib)) {
    if (ParSib %in% c("par", "sib", "ped")) {
      Module <- switch(ParSib, par = "par", sib = "ped", 
                       ped = "ped")
    }
    else {
      stop("'ParSib' must be 'par' or 'sib'; or use only 'Module' instead")
    }
  }
  if (!Module %in% c("par", "ped")) 
    stop("'Module' must be 'par' or 'ped'")
  if (!quiet) 
    message("Searching for non-assigned ", c(par = "parent-offspring", 
                                             ped = "relative")[Module], " pairs ...", " (Module = ", 
            Module, ")")
  if (!is.null(SeqList)) {
    if (!inherits(SeqList, "list")) 
      stop("'SeqList' must be a list")
    NewName = c(Pedigree = "Pedigree", PedigreePar = "Pedigree", 
                LifeHist = "LifeHistData", AgePriors = "AgePrior")
    for (x in names(NewName)) {
      if (x %in% names(SeqList)) {
        if (x == "PedigreePar" & "Pedigree" %in% names(SeqList)) 
          next
        if (!quiet) 
          message("using ", x, " in SeqList")
        assign(NewName[x], SeqList[[x]])
      }
    }
  }
  if (!quiet %in% c(TRUE, FALSE)) 
    stop("'quiet' must be TRUE or FALSE")
  GenoM <- sequoia:::CheckGeno(GenoM, quiet = quiet, Plot = FALSE)
  gID <- rownames(GenoM)
  Pedigree <- sequoia:::PedPolish(Pedigree, gID, DropNonSNPd = FALSE, 
                                  NullOK = TRUE)
  LH <- sequoia:::CheckLH(LifeHistData, gID, sorted = TRUE)
  if (!is.null(AgePrior)) {
    AgePrior <- sequoia:::CheckAP(AgePrior)
  }
  else {
    AgePrior <- sequoia:::MakeAgePrior(Pedigree, LifeHistData, Plot = FALSE, 
                                       quiet = quiet)
  }
  if ("Specs" %in% names(SeqList)) {
    if (!quiet) 
      message("settings in SeqList$Specs will overrule input parameters")
    SeqList$Specs$Module <- Module
    PARAM <- sequoia:::SpecsToParam(SeqList$Specs, SeqList$ErrM, ErrFlavour, 
                                    dimGeno = dim(GenoM), Module, MaxPairs, quiet)
  }
  else {
    PARAM <- sequoia:::namedlist(dimGeno = dim(GenoM), MaxPairs, Err, 
                                 ErrFlavour, Tfilter, Tassign, nAgeClasses = nrow(AgePrior), 
                                 MaxSibshipSize = max(table(Pedigree$dam), table(Pedigree$sire), 
                                 90, na.rm = TRUE) + 10, Module = as.character(Module), 
                                 Complex, Herm, quiet)
    PARAM$ErrM <- sequoia:::ErrToM(Err, flavour = ErrFlavour, Return = "matrix")
  }
  if (!"MaxMismatchV" %in% names(PARAM)) {
    sts <- sequoia:::SnpStats(GenoM, Plot = FALSE)
    PARAM$MaxMismatchV <- setNames(sequoia:::CalcMaxMismatch(Err = PARAM$ErrM, 
                                                             MAF = sts[, "AF"], ErrFlavour = PARAM$ErrFlavour, 
                                                             qntl = 0.1^(1/nrow(GenoM))), c("DUP", "OH", "ME"))
  }
  if (!is.na(MaxMismatch)) {
    PARAM$MaxMismatchV["OH"] <- MaxMismatch
    PARAM$MaxMismatchV["ME"] <- MaxMismatch
  }
  sequoia:::CheckParams(PARAM)
  if (any(LifeHistData$Sex == 4) && PARAM$Herm == "no") {
    if (!quiet) 
      message("detected hermaphrodites (sex=4), changing Herm to 'A'")
    PARAM$Herm <- "A"
  }
  FortPARAM <- sequoia:::MkFortParams(PARAM, fun = "mayberel")
  PedN <- sequoia:::PedToNum(Pedigree, gID, DoDummies = "new")
  LHF <- sequoia:::orderLH(LH, gID)
  
  TMP <- .Fortran(sequoia:::findambig, 
                  ng = as.integer(nrow(GenoM)), 
                  specsint = as.integer(FortPARAM$SpecsInt), 
                  specsintamb = as.integer(FortPARAM$SpecsIntAmb), 
                  specsdbl = as.double(FortPARAM$SpecsDbl), 
                  errv = as.double(FortPARAM$ErrM), 
                  genofr = as.integer(GenoM), 
                  sexrf = as.integer(LHF$Sex), 
                  byrf = as.integer(c(LHF$BirthYear, LHF$BY.min, LHF$BY.max)), 
                  aprf = as.double(AgePrior), 
                  parentsrf = as.integer(PedN$PedPar), 
                  dumparrf = as.integer(PedN$DumPar), 
                  namb = as.integer(0), 
                  ambigid = integer(2 * MaxPairs), 
                  ambigrel = integer(2 * MaxPairs), 
                  ambiglr = double(2 * MaxPairs), 
                  ambigoh = integer(MaxPairs), 
                  ntrio = as.integer(0), 
                  trioids = integer(3 * nrow(GenoM)), 
                  triolr = double(3 * nrow(GenoM)), 
                  triooh = integer(3 * nrow(GenoM)))
  
  TMP$ambiglr[abs(TMP$ambiglr - 999) < 0.1] <- NA
  TMP$ambiglr <- round(TMP$ambiglr, 2)
  TMP$ambigoh[TMP$ambigoh < 0] <- NA
  TMP$triolr[abs(TMP$triolr - 999) < 0.1] <- NA
  TMP$triolr <- round(TMP$triolr, 2)
  TMP$triooh[TMP$triooh < 0] <- NA
  
  if (TMP$namb > 0) {
    RelName <- c("PO", "FS", "HS", "GP", "FA", "HA", "U ", "Q", "2nd")
    Na <- TMP$namb
    TMP$ambigid <- sequoia:::NumToID(TMP$ambigid, 0, gID, NULL)
    AmbigRel <- factor(TMP$ambigrel, levels = 1:9, labels = RelName)
    MaybeRel <- data.frame(sequoia:::VtoM(TMP$ambigid, Na), 
                           sequoia:::VtoM(AmbigRel, Na), 
                           sequoia:::VtoM(TMP$ambiglr, Na), 
                           stringsAsFactors = FALSE)
    names(MaybeRel) <- c("ID1", "ID2", "Relx", "TopRel", "LLR_Rx_U", "LLR")
    MaybeRel <- MaybeRel[, -which(names(MaybeRel) %in% c("Relx", "LLR_Rx_U"))]
    MaybeRel$OH <- TMP$ambigoh[seq_len(Na)]
    LH$BirthYear[LH$BirthYear < 0] <- NA
    MaybeRel <- merge(MaybeRel, setNames(LH[, 1:3], c("ID1", "Sex1", "BirthYear1")), all.x = TRUE)
    MaybeRel <- merge(MaybeRel, setNames(LH[, 1:3], c("ID2", "Sex2", "BirthYear2")), all.x = TRUE)
    MaybeRel$AgeDif <- with(MaybeRel, BirthYear1 - BirthYear2)
    MaybeRel <- MaybeRel[, c("ID1", "ID2", "TopRel", "LLR", "OH", "BirthYear1", "BirthYear2", "AgeDif", "Sex1", "Sex2")]
    for (i in 1:Na) {
      if (is.na(MaybeRel$AgeDif[i])) next
      if (MaybeRel$AgeDif[i] < 0) {
        tmpRel <- MaybeRel[i, ]
        tmpRel$AgeDif <- abs(tmpRel$AgeDif)
        MaybeRel[i, ] <- tmpRel[, c("ID2", "ID1", "TopRel", "LLR", "OH", "BirthYear2", "BirthYear1", "AgeDif", "Sex2", "Sex1")]
      }
    }
    if (Module == "ped") {
      MaybeRel <- with(MaybeRel, MaybeRel[TopRel %in% c("PO", "FS", "HS", "GP", "FA", "2nd", "Q"), ])
    }
    MaybeRel <- MaybeRel[order(ordered(MaybeRel$TopRel, levels = RelName), -MaybeRel$LLR), ]
    CalcSnpdBoth <- function(Pairs, GenoM) {
      sapply(seq_along(Pairs[, 1]), function(i, G = GenoM) {
        sum(G[Pairs[i, 1], ] >= 0 & G[Pairs[i, 2], ] >= 0)
      })
    }
    if (nrow(MaybeRel) == 0) {
      MaybeRel <- NULL
    } else {
      rownames(MaybeRel) <- 1:nrow(MaybeRel)
      MaybeRel$SNPdBoth <- CalcSnpdBoth(MaybeRel[, c("ID1", "ID2")], GenoM)
    }
  } else {
    MaybeRel <- NULL
  }
  
  if (quiet < 1) {
    if (!is.null(MaybeRel)) {
      nRel <- nrow(MaybeRel)
      nPO <- sum(MaybeRel$TopRel == "PO")
    } else {
      nRel <- 0
      nPO <- 0
    }
    message("Found ", nPO, " likely parent-offspring pairs, and ", nRel - nPO, " other non-assigned pairs of possible relatives")
  }
  
  if (TMP$ntrio > 0) {
    trios <- data.frame(sequoia:::VtoM(TMP$trioids, nr = TMP$ntrio, nc = 3), 
                        sequoia:::VtoM(TMP$triolr, nr = TMP$ntrio, nc = 3), 
                        sequoia:::VtoM(TMP$triooh, nr = TMP$ntrio, nc = 3), 
                        stringsAsFactors = FALSE)
    names(trios) <- c("id", "parent1", "parent2", "LLRparent1", "LLRparent2", "LLRpair", "OHparent1", "OHparent2", "MEpair")
    for (k in 1:3) trios[, k] <- sequoia:::NumToID(trios[, k], k - 1, gID, NULL)
    trios$SNPd.id.parent1 <- CalcSnpdBoth(trios[, c("id", "parent1")], GenoM)
    trios$SNPd.id.parent2 <- CalcSnpdBoth(trios[, c("id", "parent2")], GenoM)
    if (quiet < 1) 
      message("Found ", nrow(trios), " parent-parent-offspring trios")
  } else {
    trios <- NULL
  }
  
  if (Module == "par") {
    return(list(MaybePar = MaybeRel, MaybeTrio = trios))
  } else {
    return(list(MaybeRel = MaybeRel, MaybeTrio = trios))
  }
}
