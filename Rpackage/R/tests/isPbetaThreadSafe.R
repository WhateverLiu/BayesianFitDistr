

Rcpp::sourceCpp("cpp/tests/isPbetaThreadSafe.cpp")
N = 1e7
a = runif(N, 1, 10)
b = runif(N, 1, 10)
q = runif(N)
myRst = checkPbetaThreadSafe(
  q, a, b, maxCore = 16, grainSize = 4)
actualRst = pbeta(q, a, b, log.p = T)
range(actualRst - myRst)


# Use q as p since it's in [0, 1]
N = 1e6
p = runif(N, 0.1, 0.9)
a = runif(N, 1, 10)
b = runif(N, 1, 10)
myRst = checkQbetaThreadSafe(
  p, a, b, maxCore = 2, grainSize = 4)
actualRst = qbeta(p, a, b, log.p = F)
range(actualRst - myRst)
# d = actualRst - myRst
# range(d[!is.nan(d)])








