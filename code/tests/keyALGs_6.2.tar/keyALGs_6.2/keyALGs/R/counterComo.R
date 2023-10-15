


#'
#' Counter-comonotonization
#' 
#' Compute the sum of two distributions assuming minimum correlation.
#' 
#' @inheritParams como
#' 
#' @details Countermonotonicity is undefined in three or more dimensional space. 
#' The function is equivalent to reversing the row order of \code{Y} and calling 
#' \code{\link{como}(X, Y)}.
#' 
#' @example inst/examples/counterComo.R
#'
counterComo <- function (X, Y, N = 64L, rgMethod = "lr")
{
  suppressMessages(require(data.table))
  rst = como(X, Y[nrow(Y):1, ], N, "", TRUE)
  colnames(rst) = c("val", "P")
  data.table::setDT(rst)
  rst = rst[, .(P = sum(P)), by = val]
  data.table::setDF(rst)
  rst = rst[order(rst$val), ]
  rownames(rst) = NULL
  if (N <= 1 || rgMethod == "") return(rst)
  ngd = seq(rst$val[1], rst$val[nrow(rst)], len = N)
  if (rgMethod == "lr") return(rglr(rst, ngd))
  rglr4split(rst, ngd)
}





