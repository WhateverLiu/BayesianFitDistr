





rm(list = ls())


tmp = list.files('R')
tmp = setdiff(tmp[tools::file_ext(tmp) == 'R'], 'RcppExports.R')
tmp = paste0('R/', tmp)
for (x in tmp) source(x)


tmp = list.files('src')
tmp = setdiff(tmp[tools::file_ext(tmp) == 'cpp'], 'RcppExports.cpp')
tmp = paste0('src/', tmp)
for (x in tmp) Rcpp::sourceCpp(x, cacheDir = "../tempFiles")


oldDistTable = data.table::setDF(data.table::fread(
  "../tempFiles/CovA_nPts64-debora.csv")) # For comparing ignorance score.
load("../data/cvgAdata.Rdata") # Coverage A
# load("../data/cvgCdata.Rdata") # Coverage C


# A claim’s damage ratio is considered to be zero if it is <= 0.002


oldDistTable = t(oldDistTable)
dimnames(oldDistTable)= NULL


dat = cvgData
dat$CDR[dat$CDR <= 0.002] = 0
datResv = dat[order(dat$MDR, dat$CDR), , drop = F]
rm(cvgData); gc()




if (F)
{
  
  
  tmp2 = new.env()
  load("../tempFiles/tmp.Rdata", envir = tmp2)
  tmp3 = LBFGSBtrbFitList(abc = tmp2$abc, lm1 = tmp2$lm1, 
                          empDistrList = tmp2$condEmpDistrs, 
                          abcLB = tmp2$abcLB, abcUB = tmp2$abcUB,
                          maxCore = 100,
                          sequentialUpdate = 5820
                          # sequentialUpdate = -1
  )
  
  
  
}




# Run full fit.
if (T)
{
  
  
  targetMDRs = c(1:10000, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5
  
  
  optRst = fullFit(
    mdr = dat$MDR,
    cdr = dat$CDR,
    windowSize = as.integer(round(nrow(dat) * 0.01)),
    slidingSpeed = as.integer(round(nrow(dat) * 1e-4)),
    sampleSize = round(nrow(dat) * (1 - 1 / exp(1))),
    NsampleSets = 100L,
    interpolationMethod = "linear",
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
    abc = matrix(c(4, 1, 2)),
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
  d = solve_d(param, eps= 1e-8, maxit = 100)
  abcd = param
  abcd[4, ] = d
  
  
  maxQuantile = 0.99
  maxLowerBound = 0.05
  tailProbThreshold = 1 - maxQuantile
  maxes = pmin(1, actuar_qtrbeta(
    1 - tailProbThreshold / (1 - optRst$P0), abcd))
  maxes[!is.finite(maxes)] = 0
  maxes = pmax(maxLowerBound, maxes)
  ind = longestNonDecreasingSubseq(maxes)
  maxFun = splinefun(x =  targetMDRs[ind], y = maxes[ind], method = "hyman")
  maxes = maxFun(targetMDRs)
  
  
  TrBtable = rbind(targetMDRs, maxes, optRst$P0, abcd)
  
  
  pmftable = NGFMfitDistr::fullDiscretize(
    TrBtable = TrBtable, 
    supportSizes = c(42L, 64L), 
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
  
  
  tb = as.data.frame(pmftable$tuned$p42)
  colnames(tb) = c("MDR", "max", paste0("P", 1:(ncol(tb) - 2L)))
  data.table::fwrite(tb, row.names = TRUE, 
                     file = "../tempFiles/examplePMF42.csv")
  
  
  tb = as.data.frame(pmftable$tuned$p64)
  colnames(tb) = c("MDR", "max", paste0("P", 1:(ncol(tb) - 2L)))
  data.table::fwrite(tb, row.names = TRUE, 
                     file = "../tempFiles/examplePMF64.csv")
  
  
  tb = as.data.frame(t(TrBtable))
  colnames(tb) = c("MDR", "max", "P0", "a", "b", "c", "d")
  data.table::fwrite(tb, file = "../tempFiles/exampleTrBtable.csv")
  
  
  
}




if (F)
{
  Rcpp::sourceCpp("../code/src/amalgamation.cpp", verbose = T, 
                  cacheDir = "../tempFiles")
  
  
  optRst = LBFGSBtrbFitList(
    abc = tmp$abc, 
    lm1 = tmp$lm1,
    empDistrList = tmp$empDistrList,
    abcLB = tmp$abcLB, abcUB = tmp$abcUB, 
    scaleEps = 1e-8, scaleMaxit = 100, 
    distanceFun = "likelihood",
    max_iterations = 100, maxCore = 1, 
    RIBlib = "Numerical Recipes", 
    sequentialUpdate = tmp$sequentialUpdate, 
    hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
    epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
    max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
    ftol = 1e-4, wolfe = 0.9)
  
  
}



srcFiles = list.files('src')
srcFiles = srcFiles[tools::file_ext(srcFiles) == "cpp" & 
                      srcFiles != "RcppExports.cpp"]
srcFiles = paste0("src/", srcFiles)
tmp = sapply(srcFiles, function(x) 
  Rcpp::sourceCpp(x, verbose = 1, cacheDir = '../tempFiles'))


rfiles = list.files('R', full.names = T)
rfiles = rfiles[rfiles != "R/RcppExports.R"]
rfiles = rfiles[tools::file_ext(rfiles) == "R"]
tmp = sapply(rfiles, function(x) source(x))





