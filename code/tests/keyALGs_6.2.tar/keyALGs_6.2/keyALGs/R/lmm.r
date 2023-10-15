
#'
#' Regriding via local (up to 2nd) moment matching
#' 
#' Regrid a PMF via local moment matching (LMM). The algorithm usually but not 
#' always preserves the PMF's second moment.
#'
#' @inherit rglr
#'
#' @details See \href{https://www.mdpi.com/2227-9091/7/2/54/htm}{Direct and 
#' Hierarchical Models for Aggregating Spatially Dependent Catastrophe Risks}.
#' \code{\link{LMMdiscreteToDiscrete}} would more likely inflate variance than 
#' four-point regriding \code{\link{rglr4split}}.
#' 
#' @example inst/examples/LMMdiscreteToDiscrete.R
LMMdiscreteToDiscrete <- function(X, ngd)
{
  lmm2ndM(X, ngd)
}


#' Discretize probability density function (PDF)
#' 
#' Discretize PDF using local (up to 2nd) moment matching
#' 
#' @param pdfun Probability density function (PDF) \code{pdfun(x, ...)},
#' e.g., \code{dnorm(x, ...)}.
#'
#' @param \dots Additional parameters to \code{pdfun()}.
#' 
#' @param gd A sorted numeric vector. Support of the result PMF.
#' 
#' @param subdivisions The range of \code{gd} is discretized into 
#' \code{subdivisions} intervals for numeric integration.
#'
#' @details Let \code{Min = gd[1], Max = gd[length(gd)], d = (Max - Min) / subdivisions}.
#' This function (i) evaluates the probability densities at support points 
#' \code{{Min, Min + d, Min + 2d, ..., Max}} using \code{pdfun}, (ii) creates a 
#' PMF of the support and normalized densities, (iii) regrid the PMF to 
#' \code{gd} using local moment matching. 
#' 
#' If the CDF (\code{cdfun})
#' corresponding to \code{pdfun} is computationally efficient, e.g., 
#' \code{pnorm} corresponding to \code{dnorm}, one can also (i) evaluate the cumulative 
#' probabilities at support points \code{{Min, Min + d, Min + 2d, ..., Max}} 
#' using \code{cdfun}, (ii) creates a PMF of the support and the differenced 
#' cumulative probabilities, (iii) regrid the PMF to \code{gd} using local 
#' moment matching.
#' 
#' @inherit rglr return
#' 
#' @note One should not expect the result PMF to have precisely equal moments 
#' to the theoretical distribution's, especially if the theoretical distribution 
#' has infinite support. See Examples.
#' 
#' @example inst/examples/LMMcontinousToDiscrete.R
LMMcontinousToDiscrete <- function(pdfun, ..., gd, subdivisions = 1e6L)
{
  sp = seq(gd[1], gd[length(gd)], len = subdivisions + 1L)
  tmp = data.frame(val = sp, P = pdfun(sp, ...))
  tmp$P = tmp$P / sum(tmp$P)
  lmm2ndM(tmp, gd)
}


























