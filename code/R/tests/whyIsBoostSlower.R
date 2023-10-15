

Rcpp::sourceCpp("cpp/tests/whyIsBoostSlower.cpp")


set.seed(123)
N = 100000L
q = runif(N)
a = runif(N, 0, 10)
b = runif(N, 0, 10)


system.time({
  Rrst = RIBF(q, a, b, useboost = F)
})


system.time({
  boostRst = RIBF(q, a, b, useboost = T)
})


range(Rrst - boostRst)





















