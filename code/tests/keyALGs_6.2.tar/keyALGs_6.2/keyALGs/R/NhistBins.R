

#' Optimal number of bins in a histogram
#' 
#' Compute the optimal number of bins in a histogram for
#' summarizing a sample based on the Freedman-Diaconis rule. 
#' 
#' @param x  A vector as the sample.
#'
#' @return  An integer.
NhistBins <- function(x)
{
  h = diff(quantile(x, c(0.25, 0.75))) * 2 * length(x) ^ (-1 / 3)
  as.integer(round(diff(range(x)) / h))
}













