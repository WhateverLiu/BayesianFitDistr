p = runif(100)
abcd = matrix(runif(4 * 100, 1, 10), ncol = 100)
NGFMfitDistr::actuar_qtrbeta(p, abcd)