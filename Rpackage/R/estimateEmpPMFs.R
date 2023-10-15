

#' Compute empirical PMFs
#' 
#' Estimate the empirical PMFs associated to the MDRs of interest.
#' 
#' 
#' @inheritParams ensembleP0
#'
#' @param windows  A list of size-2 integer vectors. Most likely a result from
#' from function \code{\link{computeWindows}}.
#' 
#' @param P0  A numeric vector. Modeled P0s corresponding to \code{targetMDRs}. 
#' 
#' @param regridMethod  Regrid method. \code{"lr"} implies linear regrid.
#' \code{"lmm"} implies local moment matching. \code{"r4"} implies four-point 
#' regrid. Default to \code{"lr"}.
#' 
#' @param nonDegenerateOldDistrs  A list of old PMFs for Bayesian update. The list
#' should have the same length as \code{targetMDRs}. Each PMF is a list or 
#' data frame of size 2: the 1st vector is the support, and the 2nd vector 
#' is the probabilities. Default to \code{NULL}, which implies no Bayesian update.
#' 
#' @param sampleWeightOnOldDistrs  Sampling weight in \code{[0, 1]} on the old 
#' distributions. Ignored if \code{nonDegenerateOldDistrs} is \code{NULL}.
#' 
#' @param empDistrSupportSize  Support size in each of the result PMFs.
#' 
#' @param maxCore  The number of CPU cores to be used.
#' 
#' @details More details can be found in 
#' \code{vignette("slides", package = "NGFMfitDistr")}.
#' 
#' @return A list of two objects:
#' \itemize{
#' \item{\code{biased}}{  List of empirical PMFs.}
#' \item{\code{biasedCorrected}}{  List of empirical PMFs whose supports have been
#' scaled such that their means match the prescribed MDRs. This list of PMFs
#' is not recommended to be fitted by TrB unless the biases are extreme.}
#' }
#' 
#' @inherit computeWindows details
#' 
#' @example inst/examples/estimateEmpPMFs.R
estimateEmpPMFs = function(
    mdr, 
    cdr, 
    windows, 
    P0,
    targetMDRs = c(1:1e4L, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5,
    regridMethod = "lr",
    nonDegenerateOldDistrs = NULL,
    sampleWeightOnOldDistrs = 0.5,
    empDistrSupportSize = 64L,
    maxCore = parallel::detectCores()
    )
{
  if (length(P0) != length(targetMDRs))
    stop("length(P0) != length(targetMDRs)")
  
  
  if (length(mdr) != length(cdr))
    stop("length(mdr) != length(cdr)")
  
  
  if ( min(diff(mdr)) < 0 ) stop(paste0(
    "mdr has not been sorted. ",
    "Note that cdr should also have been ordered by mdr."))
  
  
  oldDistrs = nonDegenerateOldDistrs
  if (!is.null(oldDistrs))
  {
    if (length(oldDistrs) != length(targetMDRs))
      stop(paste0(
        "List of old distributions is given, ", 
        "but the size of the list is unequal to that of targetMDRs."))
    oldDistrsMeans = unlist(lapply(oldDistrs, function(x) sum(x[[1]] * x[[2]])))
    if (max(abs(oldDistrsMeans - targetMDRs)) > 1e-6)
      stop(paste0("max(abs(oldDistrsMeans - targetMDRs)) > 1e-6. ", 
                  "The old distributions' means should have negligible ",
                  "differences from targetMDRs."))
  }
  
  
  windows = matrix(unlist(windows), nrow = 2)
  
  
  empDistrs = makeEmpDistrList( # Without P0.
    X = cdr, windows = windows, rstSizes = 255L, 
    regridMethod = regridMethod, maxCore = maxCore, 
    fixedMin = 1e300, fixedMax = 1e300)
  
  
  m = apply(windows, 2, function(x) mean(mdr[x[1]:x[2]]))
  if (m[1] <= 0) 
    stop("Window has nonpositive mean MDR. Please remove this window of data.")
  
  
  rawP0model = NULL
  targetEmpDistrsP0 = P0
  
  
  mainTargetMDRs = targetMDRs[
    targetMDRs >= m[1] & targetMDRs <= m[length(m)]]
  targetEmpDistrsMain = findEmpDistrGivenMDR( # No degenerate PMFs.
    empDistrs, m, mainTargetMDRs, 255L, maxCore, regridMethod) # mainTargetMDRs is MDRwanted.
  
  
  # ============================================================================
  # Use linear scaling for extrapolating empirical PMFs whose target MDRs are
  #   lower than the minimum MDR in the given data.
  # ============================================================================
  lowerInds = which(targetMDRs < m[1])
  tmpDistr = targetEmpDistrsMain[[1]]
  r = targetMDRs[lowerInds] / targetMDRs[lowerInds[length(lowerInds)] + 1L]
  lowerDistrs = lapply(r, function(mpr)
  {
    list(val = tmpDistr[[1]] * mpr, P = tmpDistr[[2]] + 0.0)
  })
  
  
  # ============================================================================
  # Use linear scaling for extrapolating empirical PMFs whose target MDRs are
  #   greater than the maximum MDR in the given data.
  # ============================================================================
  upperInds = which(targetMDRs > m[length(m)])
  tmpDistr = targetEmpDistrsMain[[length(targetEmpDistrsMain)]]
  r = targetMDRs[upperInds] / targetMDRs[upperInds[1] - 1L]
  upperDistrs = upScaleAndBoundPMF(
    pmf = tmpDistr, upperBound = 1.0, upscaler = r, regridMethod = regridMethod)
  
  
  targetEmpDistrs = c(lowerDistrs, targetEmpDistrsMain, upperDistrs)
  
  
  if (length(targetEmpDistrs) != length(P0))
    stop("length(targetEmpDistrs) != length(P0)")
  
  
  includeP0 = function(distlist)
  {
    for (i in 1:length(distlist))
    {
      distlist[[i]][[1]] = c(0, distlist[[i]][[1]])
      distlist[[i]][[2]] = c(P0[i], distlist[[i]][[2]] * (1 - P0[i]))
    }
    distlist
  }
  
  
  if (is.null(oldDistrs)) # No old distribution for updating.
  {
    targetEmpDistrs = mixDistrList(
      targetEmpDistrs, targetEmpDistrs,
      empDistrSupportSize - 1L, 0, regridMethod, maxCore)
    
    
    mdr2impose = targetMDRs / (1.0 - P0)
    targetEmpDistrsCorrected = correctMeanBias(
      distlist = targetEmpDistrs, mdrs = mdr2impose, lim = 1.0,
      eps = 1e-10, maxCore = maxCore)
    
    
    return(list(biased = includeP0(targetEmpDistrs), 
                biasCorrected = includeP0(
                  targetEmpDistrsCorrected)))
  }
  
  
  oldDistrsP0 = unlist(lapply(oldDistrs, function(x) x[[2]][1]))
  w = sampleWeightOnOldDistrs
  P0 = (1 - w) * P0 + w * oldDistrsP0
  oldDistrWithoutP0 = lapply(oldDistrs, function(x)
  {
    list(x[[1]][-1], x[[2]][-1] / (1 - x[[2]][1]))
  })
  
  
  targetEmpDistrs = mixDistrList(
    targetEmpDistrs, oldDistrWithoutP0, empDistrSupportSize - 1L,
    sampleWeightOnOldDistrs, regridMethod, maxCore)
  
  
  mdr2impose = targetMDRs / (1.0 - P0)
  targetEmpDistrsCorrected = correctMeanBias(
    distlist = targetEmpDistrs, mdrs = mdr2impose, lim = 1.0,
    eps = 1e-10, maxCore = maxCore)
  
  
  list(biased = includeP0(targetEmpDistrs), 
       biasCorrected = includeP0(targetEmpDistrsCorrected))
}















