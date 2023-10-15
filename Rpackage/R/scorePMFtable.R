#'
#' Score PMF Table
#' 
#' Given claims data, assess the goodness of fit of the PMF table by computing
#' (i) ignorance scores and (ii) negative log-likelihood of the claims data.
#' 
#' @param pmfTable  PMF table as a numeric matrix. Every column is a 
#' (MDR, max, P0, P1, ..., Pmax).
#' 
#' @param claimsData  Claims data as a 2-column data frame. The first column
#' is MDR, and the second is claim damage ratio. The MDR should be 
#' \eqn{\in(0, 1)}. Claim damage ratio should be \eqn{\in[0, 1]}.
#' 
#' @details For computing the negative log-likelihood of the claims data, the
#' PMF table is seen as a 2D joint PMF, and claims data (MDR, claim damage ratio)
#' are realizations of the joint. More details can be found in
#' \href{../doc/slides.pdf}{\code{vignette("slides", package = "NGFMfitDistr")}}
#'
#' @return A list of 2 objects:
#' \itemize{
#'   \item{\code{pmfsAndScores}}{ 
#'     A list. The size is the number of PMFs that have been scored. 
#'     \code{pmfsAndScores[[i]]} is a list of 3 objects:
#'     \itemize{
#'       \item{\code{scoredPMF}}{ 
#'         The scored PMF as a list of 2 numeric vectors: 
#'         support and probability.
#'       }
#'       \item{\code{claims}}{ The claim damage ratios used to score this PMF.}
#'       \item{\code{ignoranceScores}}{ Ignorance score from each damage ratio.}
#'     }
#'   }
#'   \item{\code{nllh}}{  Negative log-likelihood of the claims data. }
#' }
scorePMFtable = function(pmfTable, claimsData)
{
  
  
  dat = claimsData
  colnames(dat) = c("MDR", "CDR")
  
  
  if (!is.finite(dat$MDR) | is.finite(dat$CDR))
    stop("Non finite values exist in claimsData for scoring.")
  
  
  tmp = dat$MDR <  0.1; dat$MDR[tmp] = round(dat$MDR[tmp], 5)
  tmp = dat$MDR >= 0.1; dat$MDR[tmp] = round(dat$MDR[tmp], 4)
  
  
  if (sum( dat$MDR >= 1 | dat$MDR <= 0 | dat$CDR > 1 | dat$CDR < 0)) 
    stop(paste0(
      "Rounded MDRs are not in (0, 1) or claim damage ratios are not in [0, 1]. ",
      "Note that MDR below 0.1 are rounded to 5 decimals, otherwise 4 decimals."))
  
  
  X = pmfTable
  tableMDRint = as.integer(round(X[[1]][1, , drop = F] * 1e5))
  
  
  lisOfClaims = as.list(aggregate(
    list(CDR = dat$CDR), 
    list(ind = as.integer(round(dat$MDR * 1e5))), 
    function(x) x, simplify = FALSE))
  lisOfClaims$whichPMFs = match(lisOfClaims$ind, tableMDRint)
  lisOfClaims$ind = NULL
  
  
  pmfsBeingScored = X[, lisOfClaims$whichPMFs, drop = F]
  pmfsBeingScored = apply(pmfsBeingScored, 2, function(x)
  {
    P = x[-c(1, 2)]
    list(val = seq(0, x[2], len = length(P)), P = P)
  })
  scores = mapply(function(x, y)
  {
    scoredPMF = x
    claims = y
    s = igscore(scoredPMF, y, redistributePwithin = T)
    list(scoredPMF = scoredPMF, 
         claims = claims, 
         ignoranceScores = s)
  }, pmfsBeingScored, lisOfClaims$CDR, SIMPLIFY = F)
  
  
  sumScoreOfEachPMF  = unlist(lapply(scores, function(x) sum(x$ignoranceScore)))
  nllh = sum(sumScoreOfEachPMF) + nrow(dat) * log2(ncol(pmfTable) - 2L)
  nllh = nllh * log(2)
  
  
  list(pmfsAndScores = scores,
       nllh = nllh)
}



















