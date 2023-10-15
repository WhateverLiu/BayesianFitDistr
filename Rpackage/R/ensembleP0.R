modelP0 = function(mdr, cdr, windows, targetMDRs, 
                   eps = 1e-10, P0modeltype = "hyman", 
                   returnInfo = FALSE, warnUsers = TRUE)
{
  rst = matrix(unlist(lapply(windows, function(x)
  {
    ind = x[1]:x[2]
    r = cdr[ind]
    c(mean(mdr[ind]), sum(r < eps) / length(r))
  })), nrow = 2)
  m = c(0, rst[1,], 1)
  p0 = c(1, rst[2,], 0)
  ind = longestNonDecreasingSubseq(-p0)
  m = m[ind]
  p0 = p0[ind]
  
  
  if (P0modeltype == "hyman")
    cbfun = splinefun(x = m, y = p0, method = 'hyman')
  else if (P0modeltype == "linear")
    cbfun = approxfun(x = m, y = p0)
  else stop("Interpolation method unknown.")
  P0 = cbfun(targetMDRs)
  if (warnUsers && !all(P0 < 1 - targetMDRs))
    warning(paste0(
      "P0s are too high: P0 >= 1 - targetMDR. ",
      "Only Bernoulli distributions could be fitted to the data. ",
      "Remove some zero claims to make Transformed Beta modeling feasible. ",
      "Note that this warning is generated only for one ensemble member. ",
      "The ensemble average can still be valid."))
  if (returnInfo)
  {
    rst = as.data.frame(rst)
    colnames(rst) = c("MDR", "P0")
    list(trainingData = rst, cbfun = cbfun, P0 = P0)  
  }
  else P0
}


#' 
#' Estimate \eqn{P_0}
#' 
#' Model \eqn{P_0} from data via the ensemble approach. The ensemble size can
#' be set to 1.
#' 
#' @param mdr  MDRs in claims data. \code{mdr} and \code{cdr} will be reshuffled and 
#' reordered internally.
#' 
#' @param cdr  Claim damage ratios in claims data.
#' 
#' @param windowSize  A positive integer as the sliding window size.
#' 
#' @param slidingSpeed  A positive integer as the sliding speed.
#' 
#' @param NsampleSets  Ensemble size. Default to 100.
#' 
#' @param sampleSize  Sample size for each member in the ensemble.
#' 
#' 
#' @param targetMDRs  Prescribed MDRs. Default to the 18999 MDRs in coverage
#' A, B, C PMF tables: 0.00001, 0.00002, ..., 0.09999, 0.1000, 0.1001, 
#' 0.1002 ..., 0.9999. 
#' 
#' @param interpolationMethod  Interpolation method, either \code{"linear"} or 
#' \code{"hyman"}. Default to \code{"linear"}, which invokes 
#' \code{\link[stats]{approxfun}(method = "linear")}.
#' \code{"hyman"} invokes \code{\link[stats]{splinefun}(method = "hyman")}. 
#' The argument controls how we interpolate the longest monotonic subsequence
#' of \eqn{P_0} (or \eqn{\max}). The interpolation function is then used to
#' estimate the \eqn{P_0} (or \eqn{\max}) at any MDR.
#' 
#' @param linearIntpoThenCpp  A boolean. If \code{interpolationMethod} is 
#' \code{"linear"} and \code{linearIntpoThenCpp} is true, use the speed 
#' optimized version of routine for modeling \eqn{P_0}. 
#' 
#' @param randomSeed  Random seed. Default to \code{42}.
#' 
#' @param maxCore  Maximum number of CPU cores to use. 
#' Default to \code{parallel::detectCores()}.
#' 
#' @param tempDir  A directory for saving computing progress. Default to
#' \code{"../tempFiles/CharlieTempMP/C"}. Will be created if nonexistent.
#' 
#' @param verbose  \code{TRUE} prints computing progress. Default \code{TRUE}.
#' 
#' @return A list of the \code{NsampleSets} \eqn{P_0} models.
#' 
#' @inherit computeWindows details
#' 
#' @note  If the ensemble approach is not to be used, set \code{sampleSize} to
#' \code{length(mdr)}, and set \code{NsampleSets} to \code{1}.
#' 
#' @example inst/examples/ensembleP0.R
#'
ensembleP0 = function(
  mdr,
  cdr,
  windowSize,
  slidingSpeed,
  sampleSize,
  NsampleSets = 100L,
  targetMDRs = c(1:10000, seq(1e4L+10L, 1e5L-10L, by=10L))/1e5,
  interpolationMethod = "linear",
  linearIntpoThenCpp = TRUE,
  randomSeed = 42L,
  maxCore = parallel::detectCores(),
  tempDir = "../tempFiles/CharlieTempMP/C",
  verbose = TRUE
  )
{
  
  dat = data.frame(MDR = mdr, CDR = cdr)
  if (verbose) cat("Computing the P0 ensemble...\n")
  
  
  # The original R approach.
  if (interpolationMethod != "linear" | !linearIntpoThenCpp)
  {
    commonData = list(
      datResv = dat, windowSize = windowSize, 
      slidingSpeed = slidingSpeed, sampleSize = sampleSize, 
      targetMDRs = targetMDRs, interpolationMethod = interpolationMethod)
    
    
    f = function(seed, commonData)
    {
      datResv = commonData$datResv
      windowSize = commonData$windowSize
      slidingSpeed = commonData$slidingSpeed
      sampleSize = commonData$sampleSize
      targetMDRs = commonData$targetMDRs
      interpolationMethod = commonData$interpolationMethod
      dat = datResv
      set.seed(seed)
      dat = dat[sample(nrow(dat), sampleSize), , drop = F]
      dat = dat[order(dat$MDR), , drop = F]
      windows = computeWindows(
        sampleSize, windowSize, start = 1L, speed = slidingSpeed)
      modelP0(
        mdr = dat$MDR, cdr = dat$CDR, windows = windows, 
        targetMDRs = targetMDRs, eps = 1e-10, 
        P0modeltype = interpolationMethod, warnUsers = TRUE)
    }
    
    
    p0models = CharliePara(
      X = as.list(seq(randomSeed, len = NsampleSets, by = 1L)),
      commonData = commonData, fun = f, maxNprocess = maxCore, 
      MPtempDir = tempDir)
  }
  else # The optimized approach implemented in C++.
  {
    
    p0models = p0ensembleLDSlinear(mdr, 
                                   cdr,
                                   windowSize,
                                   slidingSpeed,
                                   1e-10,
                                   targetMDRs,
                                   NsampleSets,
                                   # sampleSizeRatio = sampleSize / length(mdr),
                                   sampleSize = sampleSize,
                                   randomSeed,
                                   maxCore)
    # p0models = p0models$p0models
  }
  
  
  # print(str(p0models))
  
  
  p0models
  
}





















