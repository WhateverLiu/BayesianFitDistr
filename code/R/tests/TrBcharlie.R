

# ==============================================================================
# Test TrB logpdf, cdf and limited mean under 1.
# ==============================================================================
Rcpp::sourceCpp('cpp/tests/TrBcharlieTest.cpp', verbose = 1)
N = 10000
n = 1000
l = 1
u = 20
dat = lapply(1:N, function(x)
{
  abcd = sort(runif(4, l, u))
  abcd[1:3] = sort(abcd[1:3])
  abcd[1:2] = sample(abcd[1:2])
  abcd[4] = sample(c(runif(1, 0, 1), runif(1, 1, 10), runif(1, 10, 100)), 1)
  list(x = runif(n), abcd = abcd)
})


system.time({actuarRst = lapply(dat, function(x)
{
  abcd = x[[2]]
  a = abcd[1]; b = abcd[2]; c = abcd[3]; d = abcd[4]
  x = x[[1]]
  list(logpdf = actuar::dtrbeta(x, a / c, c, b / c, scale = d, log = T),
       cdf = actuar::ptrbeta(x, a / c, c, b / c, scale = d),
       l1m = actuar::levtrbeta(1, a / c, c, b / c, scale = d))
})})


system.time({myRst = lapply(dat, function(x)
{
  abcd = x[[2]]
  a = abcd[1]; b = abcd[2]; c = abcd[3]; d = abcd[4]
  x = x[[1]]
  testTrBcharlie(x, abcd)
})})


err = mapply(function(x, y)
{
  c(max(2 * abs((x[[1]] - y[[1]]) / (x[[1]] + y[[1]]))),
  max(2 * abs((x[[2]] - y[[2]]) / (x[[2]] + y[[2]]))),
  max(2 * abs((x[[3]] - y[[3]]) / (x[[3]] + y[[3]]))))
}, myRst, actuarRst, SIMPLIFY = F)


errInd = which.max(unlist(lapply(err, function(x) max(x))))
# dat[[errInd]]
# myRst[[errInd]]
# actuarRst[[errInd]]


range(err[[errInd]])





# ==============================================================================
# Check TrB solving for d.
# ==============================================================================
l = 1
u = 20
N = 10000L


param = lapply(1:N, function(i)
{
  abc = sort(runif(3, l, u))
  abc[1:3] = sort(abc[1:3])
  abc[1:2] = sample(abc[1:2])
  lm1 = runif(1)
  c(abc, lm1)
})


Rcpp::sourceCpp('cpp/tests/TrBcharlieTest.cpp', verbose = F)
fittedRst = lapply(param, function(x)
{
  dfromBisec = testCharlieTrbReset(
    x[1], x[2], x[3], lm1, eps = 1e-8, maxit = 100, useNewton = F)
  dfromNewton = testCharlieTrbReset(
    x[1], x[2], x[3], lm1, eps = 1e-8, maxit = 100, useNewton = T)
  bisecRst = actuar::levtrbeta(
    1, x[1] / x[3], x[3], x[2] / x[3], scale = dfromBisec)
  newtonRst = actuar::levtrbeta(
    1, x[1] / x[3], x[3], x[2] / x[3], scale = dfromNewton)
  c(bisecRst, newtonRst, lm1)
})
max(unlist(lapply(fittedRst, function(x) diff(range(x)))))
# On average, 3.5 iterations can converge.





dfromBisec = testCharlieTrbReset(
  abcd[1], abcd[2], abcd[3], lm1, eps = 1e-8, maxit = 100, useNewton = F)
dfromNewton = testCharlieTrbReset(
  abcd[1], abcd[2], abcd[3], lm1, eps = 1e-8, maxit = 100, useNewton = T)
bisecRst = actuar::levtrbeta(
  1, abcd[1] / abcd[3], abcd[3], abcd[2] / abcd[3], scale = dfromBisec)
newtonRst = actuar::levtrbeta(
  1, abcd[1] / abcd[3], abcd[3], abcd[2] / abcd[3], scale = dfromNewton)
bisecRst; newtonRst; lm1




# ==============================================================================
# Test discretization.
# ==============================================================================
Rdevel
source("R/rfuns.R")
Rcpp::sourceCpp('cpp/tests/TrBcharlieTest.cpp', verbose = F)
l = 1
u = 20


abc = sort(runif(3, l, u))
abc[1:3] = sort(abc[1:3])
abc[1:2] = sample(abc[1:2])
lm1 = runif(1)
abc_lm1 = c(abc, lm1)
# firstPoint = 1 / (64 - 1)
lastPoint  = 1
myRst = discretize(abc_lm1, lastPoint, 
                   size = 63, eps = 1e-8, maxit = 100, 
                   lastProbLowerBound = 1e-10,
                   useNewton = F)
abcd = myRst$abcd
myRst = myRst$distr
keyALGs::Mean(myRst) - lm1


myRst = discretize(abc_lm1, lastPoint, 
                   size = 63, eps = 1e-8, maxit = 100, 
                   lastProbLowerBound = 1e-10,
                   useNewton = T)
abcd = myRst$abcd
myRst = myRst$distr
keyALGs::Mean(myRst) - lm1










































































