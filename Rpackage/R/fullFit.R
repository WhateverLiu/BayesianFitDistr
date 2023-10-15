

#'
#' End-to-end Inflated Transformed Beta Fitting Pipeline
#'
#' An amalgamation of the building blocks in the fitting methodology,
#' including \code{\link{computeWindows}}, \code{\link{ensembleP0}}, 
#' \code{\link{estimateEmpPMFs}}, 
#' \code{\link{extractMain}}, \code{\link{LBFGSBtrbFitList}}, etc. 
#' The pipeline imports claims data, exports the inflated TrB parameter table,
#' manages deductibles, and incorporates Bayesian updates.
#' Outputs from major computational stages are saved in the return value
#' for diagnostics. Diagnostic visualizations can be made during the computation.
#' Users can type \code{fullFit} in the R console and read the documented code.
#' 
#' @inheritParams ensembleP0
#' 
#' @inheritParams estimateEmpPMFs
#' 
#' @inheritParams LBFGSBtrbFitList
#' 
#' @param deductible  Either \code{NULL} or a numeric vector of size 
#' \code{length(mdr)}. If not \code{NULL}, all elements should be 
#' \eqn{\in[0,1]}. Elements are ratios of deductibles to replacement values. 
#' Default to \code{NULL}.
#' 
#' @param isDedFranchise  A boolean. \code{TRUE} implies the deductibles are
#' franchise deductibles. Default to \code{TRUE}. Ignored if 
#' \code{deductible==NULL}.
#' 
#' @param givenP0  Either \code{NULL} or a numeric vector of size 
#' \code{length(mdr)}. If not \code{NULL}, \code{givenP0} will be the model 
#' for \eqn{P_0}.
#' 
#' @param figureDir  Path to the directory for saving plots. 
#' \code{figureDir = ""} silences visualization. The directory will be created
#' if nonexistent. Default to \code{"../figure"}.
#' 
#' @param fitToBiasCorrectedPMFs  A boolean. \code{TRUE} will fit TrBs to
#' bias corrected empirical PMFs. Default to \code{FALSE}.
#' 
#' @param startingPmf  An integer default to \code{0}:
#' \itemize{
#' \item{\code{-1}:}{ Invokes the plain fitting procedure, that is, for each PMF, 
#' the function will use each of the \eqn{K} initializations in \code{abc} to 
#' train a TrB, and among the \eqn{K} optimized parameters, select the one 
#' that leads to the best objective function value.}
#' 
#' 
#' \item{\code{0}:}{ Invokes the bidirectional sequential 
#' fitting procedure. It chooses the PMF associated to the median MDR in data
#' as the starting PMF. For each of the \eqn{K} initializations in \code{abc}, 
#' the function will fit the starting PMF, and then use the trained TrB 
#' parameters as the initialization for fitting the next/prior PMF. This
#' will end up with \eqn{K} optimized parameters for each PMF. The function 
#' will then return the one that leads to the best objective function value.
#' }
#' 
#' \item{\code{>=1}:}{  Similar to \code{0} but the index of 
#' the starting PMF is fixed to \code{startingPmf}. }
#' }
#'
#' @return A list of objects procedurally generated during major computational 
#' stages. All the objects except for \code{TrBtable} are intermediaries 
#' saved for diagnostics.
#' \itemize{
#'   \item{\code{maxCore}:}{  The number of cores that were actually used. }
#'   \item{\code{mdrRangeInData}:}{  The range of MDRs in data. }
#'   \item{\code{p0models}:}{  Output from \code{\link{ensembleP0}}. }
#'   \item{\code{P0}:}{  The \eqn{P_0} model. A numeric vector of size 
#'     \code{length(targetMDRs)}. }
#'   \item{\code{windows}:}{  Output from \code{\link{computeWindows}}. }
#'   \item{\code{empDistrsLists}:}{  Output from \code{\link{estimateEmpPMFs}}. }
#'   \item{\code{condEmpDistrs}:}{  Output from \code{\link{extractMain}}. }
#'   \item{\code{startingPmf}:}{  Index of the first PMF that was fitted. }
#'   \item{\code{startingMDR}:}{  MDR associated to the starting PMF. 
#'     This only exists if \code{startingPmf >= 0}.  }
#'   \item{\code{optRst}:}{  Output from \code{\link{LBFGSBtrbFitList}}. }
#'   \item{\code{postDeductibleP0}:}{  
#'     The \eqn{P_0} model after subtracting the missing probability mass due 
#'     to deductibles. This only exists if \code{deductible} is not \code{NULL}. }
#'   \item{\code{TrBtable}:}{  Inflated Transformed Beta parameter table as a 
#'     6-column dataframe. The column names are \code{MDR,P0,a,b,c,d}. 
#'     Note that \code{max} is still missing and will be determined before
#'     discretization.}
#' }
#' 
#' 
fullFit = function(
    mdr,
    cdr,
    windowSize,
    slidingSpeed,
    sampleSize,
    NsampleSets = 100L,
    interpolationMethod = "linear",
    linearIntpoThenCpp = TRUE,
    targetMDRs = c(1:10000, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5,
    randomSeed = 123L,
    maxCore = parallel::detectCores(),
    deductible = NULL,
    isDedFranchise = TRUE,
    givenP0 = NULL,
    nonDegenerateOldDistrs = NULL,
    sampleWeightOnOldDistrs = 0.8,
    empDistrSupportSize = 64L,
    regridMethod = "lr",
    tempDir = "../tempFiles/CharlieTempMP/C",
    figureDir = "../figure",
    fitToBiasCorrectedPMFs = FALSE,
    abc = matrix(c(4, 5, 6)),
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
{
  
  result = list()
  
  
  # ============================================================================
  # Set the maximum number of CPU cores that can be used.
  # ============================================================================
  maxCore = max(min(parallel::detectCores(), maxCore), 1L)
  result$maxCore = maxCore
  
  
  # ============================================================================
  # Data preprocessing. Order by MDR, CDR + random shuffle + order by MDR
  # ============================================================================
  dat = data.frame(MDR = mdr, CDR = cdr)
  # datResv = dat[order(dat$MDR, dat$CDR), , drop = F]
  dat = dat[order(dat$MDR, dat$CDR), , drop = F]
  set.seed(randomSeed)
  dat = dat[sample(nrow(dat)), , drop = F]
  dat = dat[order(dat$MDR), , drop = F]
  
  
  mdrRangeInData = c(
    mean(dat$MDR[1:(windowSize - 1L)]), 
    mean(dat$MDR[(length(dat$MDR) - windowSize + 1L):length(dat$MDR)]))
  result$mdrRangeInData = mdrRangeInData
  
  
  # ============================================================================
  # Model P0. If P0 is given by user, do nothing.
  # ============================================================================
  if (!is.null(givenP0)) P0 = givenP0 # P0 is given.
  else
  {
    # ==========================================================================
    # To disable the ensemble estimation, `sampleSize` should have been set to 
    #   `nrow(dat)` and `NsampleSets` should have been set to 1.
    # ==========================================================================
    p0models = ensembleP0(
      mdr = dat$MDR, cdr = dat$CDR, windowSize = windowSize, 
      slidingSpeed = slidingSpeed, sampleSize = sampleSize,
      NsampleSets = NsampleSets, interpolationMethod = interpolationMethod, 
      linearIntpoThenCpp = linearIntpoThenCpp,
      targetMDRs = targetMDRs, randomSeed = randomSeed + 1007L, 
      maxCore = maxCore, tempDir = tempDir, verbose = verbose)
    result$p0models = p0models
    
    
    if (figureDir != "")
    {
      dir.create(figureDir, showWarnings = FALSE)
      figPath = paste0(figureDir, "/p0models.png")
      plotP0model(figPath, targetMDRs, p0models)
      if (verbose) cat(
        "P0 plot is saved to ", normalizePath(figPath), "\n")
    }
    P0 = rowMeans(as.data.frame(p0models))
  }
  
  
  if (verbose) cat("Checking the validity of P0 ensemble mean...\n")
  
  
  if (!all(P0 < 1 - targetMDRs))
  {
    stop(
      paste0("P0s are too high: P0 >= 1 - targetMDR. Only Bernoulli distributions ",
             "could be fitted to the data. Remove some zero claims to make ",
             "Transformed Beta modeling feasible.")
    )
  } else 
  {
    if (verbose) cat("P0 model is valid.\n")
  }
  
  
  result$P0 = P0
  
  
  # ============================================================================
  # Remove zero claims from data after P0 has been modeled.
  # ============================================================================
  mdrMedian = median(dat$MDR)
  dat = dat[dat$CDR > 1e-10, , drop = F]
  
  
  # ============================================================================
  # If deductibles exist and if they are ordinary deductibles, shift the data
  #   to the right.
  # ============================================================================
  if (!is.null(deductible) & !isDedFranchise) dat$CDR = dat$CDR + deductible
    
  
  # ============================================================================
  # Compute windows.
  # ============================================================================
  windows = computeWindows(
    xlen = sampleSize, windowSize = windowSize, start = 1L, 
    speed = slidingSpeed, returnList = TRUE)
  result$windows = windows
  
  
  # ============================================================================
  # Compute the empirical PMFs corresponding to the 18999 MDRs, aka `targetMDRs`.
  # ============================================================================
  empDistrsLists = estimateEmpPMFs(
    mdr = dat$MDR, cdr = dat$CDR, windows = windows, P0 = P0, 
    targetMDRs = targetMDRs, regridMethod = regridMethod,
    nonDegenerateOldDistrs = nonDegenerateOldDistrs, 
    sampleWeightOnOldDistrs = sampleWeightOnOldDistrs, 
    empDistrSupportSize = empDistrSupportSize, 
    maxCore = maxCore)
  result$empDistrsLists = empDistrsLists
  
  
  # ============================================================================
  # There are two versions of PMFs in `empDistrsLists`: bias corrected and 
  #   the original. Very rarely do we want to choose the bias corrected for
  #   fitting TrBs.
  # ============================================================================
  if (fitToBiasCorrectedPMFs) empDistrs = empDistrsLists$biasCorrected else
    empDistrs = empDistrsLists$biased
  
  
  # ============================================================================
  # Extract the main parts, namely PMFs without P0s in `empDistrs`. The function
  #   `extractMain()` returns two objects: the normalized PMFs without P0s, 
  #   and the target limited mean for TrBs. `lm1` is deduced from MDR and P0.
  # ============================================================================
  tmp = extractMain(
    empDistrs, targetMDRs, normalizeMainPart = TRUE, removeZeroPs = FALSE)
  condEmpDistrs = tmp$conditionalDistrs # PMFs to be fitted by TrBs.
  lm1 = tmp$lm1 # Target limited means for TrBs.
  result$condEmpDistrs = condEmpDistrs
  
  
  # ============================================================================
  # Choose which PMF in `condEmpDistrs` should be fitted first. The
  #   bi-directional algorithm will fit the next/prior PMFs using the
  #   TrB parameters optimized for the current PMF. 
  # ============================================================================
  if (startingPmf == 0)
    startingPmf = which.min(abs(mdrMedian - targetMDRs))
  
  
  result$startingPmf = startingPmf
  if (verbose & startingPmf >= 0)
  {
    cat(
      "Which PMF will be fitted first: ", startingPmf, ". ",
      "The PMF has target MDR = ", targetMDRs[startingPmf], ".\n", sep = "")
    result$startingMDR = targetMDRs[startingPmf]
  }
  
  
  # return( list (
  #   P0 = P0,
  #   abc = abc, lm1 = lm1, empDistrList = condEmpDistrs,
  #   abcLB = abcLB, abcUB = abcUB, scaleEps = scaleEps, scaleMaxit = scaleMaxit, 
  #   distanceFun = distanceFun, max_iterations = max_iterations, maxCore = maxCore, 
  #   RIBlib = RIBlib, sequentialUpdate = startingPmf
  # ))
  

  # ============================================================================
  # Run bi-directional sequential fitting algorithm based on L-BFGS-B. Most
  #   arguments of `LBFGSBtrbFitList()` control the backend L-BFGS-B solver.
  #   Do not change these arguments without consulting the library 
  #   documentation listed on `LBFGSBtrbFitList()`'s help page.
  # ============================================================================
  optRst = LBFGSBtrbFitList(
    abc = abc, lm1 = lm1, empDistrList = condEmpDistrs,
    abcLB = abcLB, abcUB = abcUB, scaleEps = scaleEps, scaleMaxit = scaleMaxit, 
    distanceFun = distanceFun, max_iterations = max_iterations, maxCore = maxCore, 
    RIBlib = RIBlib, sequentialUpdate = startingPmf, 
    hgrad, centralDiff, m, epsilon,
    epsilon_rel, past, delta, max_submin,
    max_linesearch, min_step, max_step, ftol, wolfe)
  result$optRst = optRst
  
  
  # ============================================================================
  # Plot optimized TrB parameters versus MDRs. If there are more than one
  #   initializations of a, b, c, we will print NA on the plot.
  # ============================================================================
  if (figureDir != "")
  {
    figPath = paste0(figureDir, "/TrBtransitionParam.png")
    plotTrBparameters(
      figPath = figPath,
      optRst = optRst,
      startingPmf = startingPmf, 
      abc = if (ncol(abc) == 1) abc else c(NA, NA, NA)
    )
    if (verbose) cat(
      "TrB parameters plot is saved as ", normalizePath(figPath), "\n")
  }
  
  
  # ============================================================================
  # Solve parameter `d`. Assemble TrB parameter set `abcd`.
  # ============================================================================
  param = optRst$param
  d = solve_d(param, eps = scaleEps, maxit = scaleMaxit)
  abcd = cbind(t(param[1:3, , drop = F]), d)
  
  
  # ============================================================================
  # If there exist deductibles, we use the currently fitted TrB to estimate
  #   the missing probability mass between 0 and deductible. This probability
  #   mass is then subtracted from P0. Next, we (i) re-impose monotonicity over 
  #   the new P0, (ii) recompute the target limited means given P0, and (iii) 
  #   adjust parameter `d` to impose the new target limited mean.
  # ============================================================================
  if (!is.null(deductible))
  {
    subtractP = actuar_ptrbeta(deductible, abcd) * (1 - P0)
    P0 = P0 - subtractP
    tmpP0 = c(1, P0, 0)
    tmpMDR = c(0, targetMDRs, 1)
    lnd = longestNonDecreasingSubseq(-tmpP0) # Subsequence indices.
    if (interpolationMethod == "hyman") cbfun = splinefun(
      x = tmpMDR, y = tmpP0, method = "hyman")
    else if (interpolationMethod == "linear") cbfun = approxfun(
      x = tmpMDR, y = tmpP0, method = "linear")
    else stop("Interpolation method unknown.")
    P0 = cbfun(targetMDRs) # Monotonic P0.
    result$postDeductibleP0 = P0
    
    
    # ==========================================================================
    # Compute the updated target limited means for TrBs.
    # ==========================================================================
    lm1 = targetMDRs / (1 - P0)
    param[4, , drop = F] = lm1
    
    
    # ==========================================================================
    # Adjust parameter `d`.
    # ==========================================================================
    d = solve_d(param, eps = scaleEps, maxit = scaleMaxit)
    abcd = cbind(t(param[1:3, , drop = F]), d)
  }
  
  
  # ============================================================================
  # Make and return TrB parameter table. Note that maxes are missing, and will
  #   be determined before discretization.
  # ============================================================================
  rst = data.frame(targetMDRs, P0, abcd)
  colnames(rst) = c("MDR", "P0", "a", "b", "c", "d")
  rownames(rst) = NULL
  result$TrBtable = rst
  result
  
  
}














