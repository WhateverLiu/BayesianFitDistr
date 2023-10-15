setwd("/finance_develop/Charlie/BayesianFitDistr/Rpackage")
data("NewZealandCvgA", package = 'NGFMfitDistr')
dat = NewZealandCvgA
str(dat)


targetMDRs = c(1:10000, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5


source("sourceAllPackgeCodes.R")


# optRst = NGFMfitDistr::fullFit(
optRst = fullFit(
  mdr = dat$MDR,
  cdr = dat$CDR,
  windowSize = as.integer(round(nrow(dat) * 0.01)),
  slidingSpeed = as.integer(round(nrow(dat) * 1e-4)),
  sampleSize = round(nrow(dat) * (1 - 1 / exp(1))),
  NsampleSets = 100L,
  interpolationMethod = "linear",
  linearIntpoThenCpp = TRUE,
  targetMDRs = targetMDRs,
  randomSeed = 42L,
  maxCore = parallel::detectCores(),
  deductible = NULL,
  isDedFranchise = TRUE,
  givenP0 = NULL,
  nonDegenerateOldDistrs = NULL,
  sampleWeightOnOldDistrs = 0.8,
  empDistrSupportSize = 64L,
  regridMethod = "lr",
  tempDir = "../tempFiles/CharlieMP/C",
  figureDir = "../tmpfigure",
  fitToBiasCorrectedPMFs = FALSE,
  abc = matrix(c(4, 1, 1.5)),
  abcLB = c(1.01, 0.1, 0.1),
  abcUB = c(30, 30, 30),
  startingPmf = 0L,
  scaleEps = 1e-8,
  scaleMaxit = 100L,
  distanceFun = "likelihood",
  max_iterations = 100L,
  RIBlib = "Numerical Recipes",
  sequentialUpdate = -1,
  hgrad = 0,
  centralDiff = TRUE,
  verbose = TRUE,
  m = 6,
  epsilon = 1e-5,
  epsilon_rel = 1e-5,
  past = 1,
  delta = 1e-10,
  max_submin = 10,
  max_linesearch = 20,
  min_step = 1e-20,
  max_step = 1e+20,
  ftol = 1e-4,
  wolfe = 0.9
)


param = optRst$optRst$param
# d = NGFMfitDistr::solve_d(
d = solve_d(
  param, eps= 1e-8, maxit = 100)
abcd = param
abcd[4, ] = d
# maxes = NGFMfitDistr::makeMax(
maxes = makeMax(
  MDR = targetMDRs, P0 = optRst$P0, abcd = abcd, tailThreshold = 1 - 0.99, 
  minMax = 0.05, interpolationMethod = "linear")


TrBtable = rbind(targetMDRs, maxes, optRst$P0, abcd)
# pmftable = NGFMfitDistr::fullDiscretize(
pmftable = fullDiscretize(
  TrBtable = TrBtable,
  supportSizes = c(42L, 64L), # Discretize TrBs on both 42 and 64-point supports.
  regridMethod = "lr",
  RIBlib = "Numerical Recipes",
  outputProbsInRows = TRUE,
  fineDiscretizationSize = 2000,
  maxCore = 1000,
  verbose = TRUE,
  downScaleSupport = FALSE,
  figureDir = "../tmpfigure",
  empDistrsLists = optRst$empDistrsLists,
  mdrRangeInData = optRst$mdrRangeInData,
  claimsDataForScoring = NULL
)






