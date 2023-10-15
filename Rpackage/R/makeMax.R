

#'
#' Make Maxes
#' 
#' Estimate proper maxes of the supports for TrB discretization.
#' 
#' @param MDR  A numeric vector of MDRs.
#' 
#' @param P0  A numeric vector of \eqn{P_0}s. Should have the same length of \code{MDR}.
#' 
#' @param abcd  A \eqn{4 \times} \code{length(MDR)} numeric matrix of TrB 
#' parameters.
#' 
#' @param tailThreshold  A tail probability for computing a tentative max.
#' This max is chosen such that the total probability of damage ratio above 
#' the max will not exceed \code{tailThreshold}.
#' 
#' @param minMax  Max for the lowest MDR. Default to 0.05.
#' 
#' @inheritParams ensembleP0
#' 
#' @details The function (i) computes a sequence of tentative maxes for 
#' the given \code{MDR}s, (ii) searches the longest nondecreasing subsequence,
#' (iii) interpolates the subsequence and estimates all the maxes for \code{MDR}.
#' 
#' @return A numeric vector of maxes corresponding to \code{MDR}.
#'
makeMax = function(MDR, P0, abcd, tailThreshold = 0.01, minMax = 0.05,
                   interpolationMethod = "linear")
{
  maxes = pmin(1, actuar_qtrbeta(1 - tailThreshold / (1 - P0), abcd))
  maxes[!is.finite(maxes)] = 0
  maxes = pmax(minMax, maxes)
  ind = longestNonDecreasingSubseq(maxes)
  if (interpolationMethod == "hyman")
    maxFun = splinefun(x = MDR[ind], y = maxes[ind], method = "hyman")
  else if (interpolationMethod == "linear")
    maxFun = approxfun(x = MDR[ind], y = maxes[ind], method = "linear")
  else stop("Unknown interpolation method.")
  maxFun(MDR)
}



















