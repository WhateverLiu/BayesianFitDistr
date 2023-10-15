
levtrbeta = function(abc_lm1, limit = 1, eps = 1e-8, maxit = 100)
{
  d = solve_d(abc_lm1, eps, maxit)
  if (!is.matrix(abc_lm1)) abc_lm1 = as.matrix(abc_lm1)
  actuar::levtrbeta(limit, abc_lm1[1, ] / abc_lm1[3, ], abc_lm1[3, ], 
                    abc_lm1[2, ] / abc_lm1[3, ], scale = d)
}


rtrbeta = function(n, abc_lm1, eps = 1e-8, maxit = 100)
{
  
  d = solve_d(abc_lm1, eps, maxit)
  if (!is.matrix(abc_lm1)) abc_lm1 = as.matrix(abc_lm1)
  actuar::rtrbeta(n, abc_lm1[1, ] / abc_lm1[3, ], abc_lm1[3, ], 
                  abc_lm1[2, ] / abc_lm1[3, ], scale = d)
}


dtrbeta = function(x, abc_lm1, eps = 1e-8, maxit = 100)
{
  d = solve_d(abc_lm1, eps, maxit)
  if (!is.matrix(abc_lm1)) abc_lm1 = as.matrix(abc_lm1)
  actuar::dtrbeta(x, abc_lm1[1,] / abc_lm1[3,], abc_lm1[3,], 
                  abc_lm1[2,] / abc_lm1[3,], scale = d)
}


ptrbeta = function(q, abc_lm1, eps = 1e-8, maxit = 100)
{
  d = solve_d(abc_lm1, eps, maxit)
  if (!is.matrix(abc_lm1)) abc_lm1 = as.matrix(abc_lm1)
  actuar::ptrbeta(q, abc_lm1[1,] / abc_lm1[3,], abc_lm1[3,], 
                  abc_lm1[2,] / abc_lm1[3,], scale = d)
}


qtrbeta = function(p, abc_lm1, eps = 1e-8, maxit = 100)
{
  d = solve_d(abc_lm1, eps, maxit)
  if (!is.matrix(abc_lm1)) abc_lm1 = as.matrix(abc_lm1)
  actuar::qtrbeta(p, abc_lm1[1,] / abc_lm1[3,], abc_lm1[3,], 
                  abc_lm1[2,] / abc_lm1[3,], scale = d)
}


actuar_rtrbeta = function(n, abcd)
{
  if (!is.matrix(abc_lm1)) abc_lm1 = as.matrix(abc_lm1)
  actuar::rtrbeta(n, abcd[1,] / abcd[3,], abcd[3,], 
                  abcd[2,] / abcd[3,], scale = abcd[4,])
}


actuar_levtrbeta = function(abcd)
{
  if (!is.matrix(abcd)) abcd = as.matrix(abcd)
  actuar::levtrbeta(1, abcd[1,] / abcd[3,], abcd[3,], 
                    abcd[2,] / abcd[3,], scale = abcd[4,])
}


actuar_dtrbeta = function(x, abcd, log = F)
{
  if (!is.matrix(abcd)) abcd = as.matrix(abcd)
  actuar::dtrbeta(x, abcd[1,] / abcd[3,], abcd[3,], 
                  abcd[2,] / abcd[3,], scale = abcd[4,], log = log)
}


actuar_ptrbeta = function(q, abcd)
{
  if (!is.matrix(abcd)) abcd = as.matrix(abcd)
  actuar::ptrbeta(q, abcd[1,] / abcd[3,], abcd[3,], 
                  abcd[2,] / abcd[3,], scale = abcd[4,])
}


#'
#' TrB Quantile
#' 
#' Inverse TrB CDF of \code{a, b, c, d}.
#'
#' @param p  A numeric vector of cumulative probabilities.
#' 
#' @param abcd  A numeric matrix of 4 rows. Each column is a vector of 
#' \code{a, b, c, d}.
#' 
#' @return A numeric vector of the quantiles, \code{q}. \code{q[i]} is the 
#' quantile corresponding to \code{p[i]} for TrB parameterized by 
#' \code{abcd[, i]}.
#'
#' @example inst/examples/actuar_qbeta.R
actuar_qtrbeta = function(p, abcd)
{
  if (!is.matrix(abcd)) abcd = as.matrix(abcd)
  actuar::qtrbeta(p, abcd[1,] / abcd[3,], abcd[3,], 
                  abcd[2,] / abcd[3,], scale = abcd[4,])
}





















