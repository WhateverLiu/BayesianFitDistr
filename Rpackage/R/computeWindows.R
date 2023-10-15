

#' Compute sliding windows
#' 
#' Compute sliding windows given window size and sliding speed.
#' 
#' @param xlen  Size of the full data.
#'
#' @param windowSize  A positive integer as the window size. Default to 1.
#' 
#' @param start  A positive integer. Index of the 1st data point included in
#' the 1st window.
#' 
#' @param speed  A positive integer. Sliding speed.
#' 
#' @param returnList  A Boolean value. If \code{TRUE}, the windows are 
#' returned as a list of integer vectors. Each integer vector is of size 2. 
#' The 1st element is the index of the first data point in window, the 2nd the
#' last data point in window. If \code{FALSE}, the list is converted to a 
#' 2-row matrix. 
#' 
#' @details More details can be found in the slide documentation 
#' \href{../doc/slides.pdf}{\code{vignette('slides', package = 'NGFMfitDistr')}}.
#' 
#' @return  If \code{returnList} is \code{TRUE}, a list of integer vectors of 
#' size 2. Otherwise an integer matrix with 2 rows.
#' 
#' 
#' @example inst/examples/computeWindows.R
computeWindows = function(
    xlen, 
    windowSize = 1L, 
    start = 1L, 
    speed = 1L, 
    returnList = TRUE
    )
{
  rst = seq(start, xlen + speed + 1L, by = speed)
  rst = rst[rst < xlen - windowSize + 1L]
  rst = c(rst, xlen - windowSize + 1L)
  if (returnList) lapply(rst, function(x) c(x, x + windowSize - 1L))
  else sapply(rst, function(x) c(x, x + windowSize - 1L))
}













