

Rcpp::sourceCpp('src/tests/testMonoInterpo.cpp', cacheDir = '../tempFiles',
                verbose = TRUE)

while(T){

  
N = 50000
# x = round(sort(runif(N)), 3)
x = sort(runif(N))
if (length(unique(x)) != N) next
# y = round(sort(runif(N)), 6)
# y = sort(runif(N))
y = round(runif(N), 5)
xtest = round(runif(N), 5)
ytest = testMonoLinearInterpo(x, y, xtest, monotype = "nonincreasing")


# tmp = NGFMfitDistr::longestIncreasingSubseq(-y)
tmp = NGFMfitDistr::longestNonDecreasingSubseq(-y)
if (length(tmp) <= 1) next
f = approxfun(x = x[tmp], y = y[tmp])
fx = f(xtest)
tmp = which(!is.na(fx))
# if (length(tmp) != 0) cat(range(fx[tmp] - ytest[tmp]))
if (length(tmp) != 0) 
{
  if (max(abs(fx[tmp] - ytest[tmp])) > 1e-5) break
}
}
cat(range(fx[tmp] - ytest[tmp]))













