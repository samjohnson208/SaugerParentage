## sequoia_custom.R by SPJ 061225
## PURPOSE: to load a custom version of sequoia::sequoia() to alter assignment/mismatch thresholds
## USAGE: Rscript sequoia_custom.R

sequoia_custom <- function (GenoM = NULL, LifeHistData = NULL, SeqList = NULL, 
                  Module = "ped", Err = 1e-04, Tfilter = -2, Tassign = 0.5, 
                  MaxSibshipSize = 100, DummyPrefix = c("F", "M"), Complex = "full", 
                  Herm = "no", UseAge = "yes", args.AP = list(Flatten = NULL, 
                  Smooth = TRUE), mtSame = NULL, CalcLLR = TRUE, quiet = FALSE, 
                  Plot = NULL, StrictGenoCheck = TRUE, ErrFlavour = "version2.0", 
                  MaxSibIter = 42, MaxMismatch = NA, FindMaybeRel = FALSE) 
{
  TimeStart <- Sys.time()
  if (!quiet %in% c(TRUE, FALSE, "verbose")) 
    stop("'quiet' must be TRUE or FALSE or 'verbose'")
  quietR <- ifelse(quiet == "verbose", FALSE, quiet)
  if (is.null(Plot)) 
    Plot <- ifelse(quietR, FALSE, TRUE)
  if (MaxSibIter != 42 && Module == "ped") {
    Module <- cut(MaxSibIter, breaks = c(-Inf, -9, -1, 0, 
                                         Inf), labels = c("pre", "dup", "par", "ped"))
    if (!quietR) 
      message("NOTE: 'MaxSibIter' will be deprecated in the future, ", 
              "please consider using 'Module' instead")
  }
  else {
    Module <- factor(Module, levels = c("pre", "dup", "par", 
                                        "ped"))
    if (is.na(Module)) 
      stop("'Module'  must be 'pre', 'dup', 'par', or 'ped'")
  }
  if (FindMaybeRel) 
    warning("'FindMaybeRel' has been deprecated and is ignored,", 
            "instead run GetMaybeRel() afterwards", immediate. = TRUE)
  if (!is.na(MaxMismatch)) 
    warning("'MaxMismatch' has been deprecated and is ignored,", 
            "now calculated automatically via CalcMaxMismatch()", 
            immediate. = TRUE)
  if (!is.null(LifeHistData) & !inherits(LifeHistData, "data.frame")) 
    stop("LifeHistData must be a data.frame or NULL")
  if (!is.null(SeqList) & !inherits(SeqList, "list")) 
    stop("SeqList must be a list or NULL")
  if (!is.null(SeqList)) {
    SeqList_names <- c("Specs", "ErrM", "args.AP", "AgePriors", 
                       "LifeHist", "PedigreePar", "MaxSibIter")
    if (!any(names(SeqList) %in% SeqList_names) | any(is.na(names(SeqList)))) 
      stop("You seem to have misspelled one or more names of elements of SeqList")
  }
  Excl <- sequoia:::CheckGeno(GenoM, quiet = quietR, Plot = Plot, Return = "excl", 
                    Strict = FALSE, DumPrefix = DummyPrefix)
  if ("ExcludedSnps" %in% names(Excl)) 
    GenoM <- GenoM[, -Excl[["ExcludedSnps"]]]
  if ("ExcludedSnps-mono" %in% names(Excl)) 
    GenoM <- GenoM[, -Excl[["ExcludedSnps-mono"]]]
  if ("ExcludedIndiv" %in% names(Excl)) 
    GenoM <- GenoM[!rownames(GenoM) %in% Excl[["ExcludedIndiv"]], 
    ]
  if ("LifeHist" %in% names(SeqList)) {
    if (!quietR) 
      message("using LifeHistData in SeqList")
    LifeHistData <- SeqList$LifeHist
  }
  else if (is.null(LifeHistData) & Module != "dup") {
    warning("no LifeHistData provided, expect lower assignment rate\n", 
            immediate. = TRUE)
  }
  ChkLH.L <- sequoia:::CheckLH(LifeHistData, gID = rownames(GenoM), sorted = FALSE, 
                     returnDups = TRUE)
  DupList <- ChkLH.L[c("DupLifeHistID", "NoLH")]
  LifeHistData <- ChkLH.L$LifeHistData
  if (!quietR) {
    gID <- rownames(GenoM)
    tbl_sex <- table(factor(LifeHistData$Sex[LifeHistData$ID %in% 
                                               gID], levels = 1:4))
    message("There are ", tbl_sex["1"], " females, ", tbl_sex["2"], 
            " males, ", tbl_sex["3"], " individuals of unkwown sex, and ", 
            tbl_sex["4"], " hermaphrodites.")
    range_Year <- matrix(NA, 4, 2, dimnames = list(c("BirthYear", 
                                                     "BY.min", "BY.max", "Year.last"), c("min", "max")))
    for (x in rownames(range_Year)) {
      range_Year[x, ] <- range(LifeHistData[, x][LifeHistData$ID %in% 
                                                   gID & LifeHistData[, x] >= 0])
    }
    message("Exact birth years are from ", range_Year[1, 
                                                      1], " to ", range_Year[1, 2])
    message("Birth year min/max are from ", min(range_Year[2:3, 
                                                           1], na.rm = TRUE), " to ", max(range_Year[2:3, 2], 
                                                                                          na.rm = TRUE))
  }
  utils::flush.console()
  if ("AgePriors" %in% names(SeqList)) {
    if (!quietR) 
      message("using AgePriors in SeqList")
    AgePriors <- sequoia:::CheckAP(SeqList[["AgePriors"]])
  }
  else {
    AgePriors <- do.call(MakeAgePrior, c(list(Pedigree = SeqList[["PedigreePar"]], 
                                              LifeHistData = LifeHistData[LifeHistData$ID %in% 
                                                                            rownames(GenoM), ], Plot = Plot, quiet = ifelse(Module == 
                                                                                                                              "dup", TRUE, quietR)), args.AP))
  }
  if ("PedigreePar" %in% names(SeqList) & Module != "dup") {
    if (!quietR) 
      message("using PedigreePar in SeqList")
    PedParents <- sequoia:::PedPolish(SeqList[["PedigreePar"]], gID = rownames(GenoM), 
                            ZeroToNA = TRUE, DropNonSNPd = TRUE)
  }
  else {
    PedParents <- NULL
  }
  utils::flush.console()
  mtDif <- sequoia:::mtSame2Dif(mtSame, gID = rownames(GenoM))
  if ("Specs" %in% names(SeqList)) {
    PARAM <- sequoia:::SpecsToParam(SeqList$Specs, SeqList$ErrM, ErrFlavour, 
                          dimGeno = dim(GenoM), Module, quiet)
    if (Module == "ped") 
      PARAM$MaxSibIter <- 42
  }
  else {
    if (is.logical(UseAge)) 
      UseAge <- ifelse(UseAge, "yes", "no")
    if ((Herm != "no" | any(LifeHistData$Sex == 4, na.rm = TRUE)) & 
        length(DummyPrefix) == 2) 
      DummyPrefix <- c(DummyPrefix, "H")
    PARAM <- sequoia:::namedlist(dimGeno = dim(GenoM), Err, ErrFlavour, 
                       Tfilter, Tassign, nAgeClasses = nrow(AgePriors), 
                       MaxSibshipSize, Module = as.character(Module), MaxSibIter, 
                       DummyPrefix, Complex, Herm, UseAge, CalcLLR, quiet)
    PARAM$ErrM <- sequoia:::ErrToM(Err, flavour = ErrFlavour, Return = "matrix")
  }
  if (!"MaxMismatchV" %in% names(PARAM)) {
    sts <- sequoia:::SnpStats(GenoM, Plot = FALSE)
    PARAM$MaxMismatchV <- setNames(sequoia:::CalcMaxMismatch(Err = PARAM$ErrM, 
                                                   MAF = sts[, "AF"], ErrFlavour = PARAM$ErrFlavour, 
                                                   qntl = 0.5^(1/nrow(GenoM))), c("DUP", "OH", "ME"))
  }
  if (any(LifeHistData$Sex == 4, na.rm = TRUE) && PARAM$Herm == 
      "no") {
    if (!quietR) 
      message("\nDetected hermaphrodites (sex=4), changing Herm to 'A'\n")
    PARAM$Herm <- "A"
  }
  if (PARAM$Herm == "B" && any(LifeHistData$Sex %in% c(1, 2))) 
    warning("Results may be inconsistent when combining Herm='B' with known-sex individuals")
  sequoia:::CheckParams(c(PARAM, list(Plot = Plot)))
  utils::flush.console()
  FortPARAM <- sequoia:::MkFortParams(PARAM, fun = "main")
  if (Module != "pre") {
    if (!quietR) 
      message("\n~~~ Duplicate check ~~~")
    DupList <- c(sequoia:::DuplicateCheck(GenoM, FortPARAM, quiet = quietR), 
                 DupList)
    utils::flush.console()
    if ("DupGenoID" %in% names(DupList)) {
      return(DupList)
    }
    else if (length(DupList) == 0 & !quietR) {
      message("No duplicates found")
    }
  }
  else DupList <- NULL
  utils::flush.console()
  if (Module == "par" | (Module == "ped" & is.null(PedParents))) {
    if (!quietR & Module == "par" & !is.null(PedParents)) {
      message("\n~~~ Parentage assignment with pedigree-prior ~~~")
    }
    else if (!quietR) {
      message("\n~~~ Parentage assignment ~~~")
    }
    ParList <- sequoia:::SeqParSib(ParSib = "par", FortPARAM, GenoM, 
                                   LhIN = LifeHistData, AgePriors = AgePriors, Parents = PedParents, 
                                   mtDif = mtDif, DumPfx = PARAM$DummyPrefix, quiet = quietR)
    if (Plot) {
      sequoia:::SummarySeq(ParList$PedigreePar, Panels = "G.parents")
    }
  }
  else if (Module != "dup" & "PedigreePar" %in% names(SeqList)) {
    ParList <- list(PedigreePar = PedParents)
  }
  else {
    ParList <- NULL
  }
  utils::flush.console()
  W <- sequoia:::tryCatch.W.E(sequoia:::getGenerations(ParList$PedigreePar, StopIfInvalid = FALSE))$warning
  if (!is.null(W)) {
    if (Module == "ped") 
      warning("Cancelling full pedigree reconstruction.")
  }
  if (Module %in% c("par", "ped") & !"AgePriors" %in% names(SeqList) & 
      is.null(W) & !"PedigreePar" %in% names(SeqList) && !is.null(LifeHistData)) {
    if (Plot) 
      Sys.sleep(1)
    AgePriors <- tryCatch(do.call(MakeAgePrior, c(list(Pedigree = ParList$PedigreePar[, 1:3], 
                                                       LifeHistData = LifeHistData[LifeHistData$ID %in% rownames(GenoM), ], 
                                                       Plot = Plot & Module == "ped", 
                                                       quiet = !(!quietR & Module == "ped")), args.AP)), 
                                                       error = function(e) {message("AgePrior error! \n", e)
                                                       return(NA)})
    if (all(is.na(AgePriors))) 
      return(ParList)
  }
  else if ("AgePriors" %in% names(SeqList) & !"PedigreePar" %in% 
           names(SeqList)) {
    if (!quietR) 
      message("using AgePriors in SeqList again")
  }
  if (nrow(AgePriors) != PARAM$nAgeClasses) {
    PARAM$nAgeClasses <- nrow(AgePriors)
    FortPARAM$SpecsInt[["nAgeCl"]] <- nrow(AgePriors)
  }
  if (Module == "ped" & is.null(W)) {
    if (!all(apply(AgePriors, 2, function(x) any(x > 0)))) 
      stop("AgePriors error: some relationships are impossible for all age differences")
    if (!quietR) 
      message("\n~~~ Full pedigree reconstruction ~~~")
    SibList <- sequoia:::SeqParSib(ParSib = "sib", FortPARAM, GenoM, 
                                   LhIN = LifeHistData, AgePriors = AgePriors, Parents = ParList$PedigreePar, 
                                   mtDif = mtDif, DumPfx = PARAM$DummyPrefix, quiet = quietR)
    ParList <- ParList[names(ParList) != "AgePriorExtra"]
    if (Plot) {
      sequoia:::PlotAgePrior(SibList$AgePriorExtra)
      sequoia:::SummarySeq(SibList$Pedigree, Panels = "G.parents")
    }
  }
  else SibList <- NULL
  if (!quietR & Module %in% c("par", "ped")) {
    message("Possibly not all ", c(par = "parents", ped = "relatives")[as.character(Module)], 
            " were assigned, consider running GetMaybeRel() conditional on this pedigree to check")
  }
  OUT <- list()
  OUT[["Specs"]] <- sequoia:::ParamToSpecs(PARAM, TimeStart, ErrFlavour)
  OUT[["ErrM"]] <- PARAM$ErrM
  if (is.function(ErrFlavour)) 
    OUT[["ErrFlavour"]] <- ErrFlavour
  if ("AgePriors" %in% names(SeqList)) {
    OUT[["args.AP"]] <- SeqList[["AgePriors"]]
  }
  else {
    OUT[["args.AP"]] <- args.AP
  }
  if (length(Excl) > 0) 
    OUT <- c(OUT, Excl)
  OUT[["AgePriors"]] <- AgePriors
  OUT[["LifeHist"]] <- LifeHistData
  if (as.numeric(Module) > 1) 
    OUT <- c(OUT, DupList)
  if (as.numeric(Module) > 2) 
    OUT <- c(OUT, ParList)
  if (as.numeric(Module) > 3) 
    OUT <- c(OUT, SibList)
  return(OUT)
}
<bytecode: 0x7ff6cc52deb8>
<environment: namespace:sequoia>