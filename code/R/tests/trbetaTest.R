

source('R/sourceCppCharlie.R')
# sourceCppCharlie('cpp/tests/trbetaTest.cpp', verbose = T)
Rcpp::sourceCpp('cpp/tests/trbetaTest.cpp', verbose = T)


N = 5L; n = 4L
x = lapply(1:N, function(x) runif(n))
abc = lapply(1:N, function(x) { x = runif(3); x[1] = x[1] + 2; x })
myRst = testTrbLimMean(x, abc)
actuarRst = mapply(function(x, abc)
{
  actuar::levtrbeta(1, abc[1] / abc[3], abc[3], abc[2] / abc[3], scale = x)
}, x, abc, SIMPLIFY = F)
range(unlist(actuarRst) - unlist(myRst))


actuarRstMean = mapply(function(x, abc)
{
  actuar::mtrbeta(1, abc[1] / abc[3], abc[3], abc[2] / abc[3], scale = x)
}, x, abc, SIMPLIFY = F)
# To let trb have finite mean, a > 1
# To let trb have finite variance, a > 2.
# To not blow the gamma function, and because a + b > 1, let c > 0.1
# The bottom line is: one needs to compute gamma((a + b) / c)
# The quadrature approx in numeric recipies is most likely unncessary for
# our case because the two parameters are vastly less than 3000.













input = lapply(1:100000, function(x) runif(4))


# Compare levtrbeta implementation.
system.time({
  actuarRst = unlist(lapply(input, function(x)
  {
    actuar::levtrbeta(1, x[1], x[2], x[3], scale = x[4])
  }))
})
system.time({
  myRst = unlist(lapply(input, function(x)
  {
    testlevtrbeta(x[1], x[2], x[3], x[4])
  }))
})
range(myRst[which(myRst < 1e100)] - actuarRst[which(myRst < 1e100)])




# Compare dtrbeta implementation
system.time({
  actuarRstDtrbeta = unlist(lapply(input, function(x)
  {
    actuar::dtrbeta(x[1], x[1], x[2], x[3], scale = x[4], log = T)
  }))
})
system.time({
  myRstDtrbeta = unlist(lapply(input, function(x)
  {
    testDtrbeta(x[1], x[1], x[2], x[3], x[4], T)
  }))
})
range(myRstDtrbeta - actuarRstDtrbeta)




# Compare ptrbeta implementation
system.time({
  actuarRst = unlist(lapply(input, function(x)
  {
    actuar::ptrbeta(x[1], x[1], x[2], x[3], scale = x[4], log.p = T)
  }))
})
system.time({
  myRst = unlist(lapply(input, function(x)
  {
    testPtrbeta(x[1], x[1], x[2], x[3], x[4], T, T)
  }))
})
range(myRst - actuarRst)








