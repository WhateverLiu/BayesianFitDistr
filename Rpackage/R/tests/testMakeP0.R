

Rcpp::sourceCpp('src/P0ensembleUsingLinearInterpo.cpp', verbose = 1,
                cacheDir = '../tempFiles')


while (T){

N = sample(1000, 1) + 1
nzero = sample(N-1, 1)
windowSize = sample(N-1, 1)
slidingSpeed = sample(N-1, 1)
eps = 1e-10
mdr = sort(runif(N))
dr = runif(N)
dr[sample(N, nzero)] = 0
rst = makeP0test(mdr, dr, windowSize, slidingSpeed, eps = eps)
rst = matrix(unlist(rst), ncol = 2)


wins = NGFMfitDistr::computeWindows(N, windowSize, start = 1L, 
                                       speed = slidingSpeed)
trueRst = t(matrix(unlist(lapply(wins, function(x)
{
  iv = x[1]:x[2]
  c(mean(mdr[iv]), sum(dr[iv] < eps) / length(iv))
})), nrow = 2))


err = max(abs(rst - trueRst))
if (err > 1e-8) break

# cat(N, nzero, windowSize, slidingspeed, "\n")
}
















