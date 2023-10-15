


Rcpp::sourceCpp("tests/pbetaTest.cpp")


set.seed(123)


N = 500000L
q = 1 - exp(runif(N, log(1e-6), log(1)))
a = exp(runif(N, log(1e-4), log(10000)))
b = exp(runif(N, log(1e-4), log(10000)))


system.time({boostRst = RIBF(q, a, b, 0)})


system.time({Rrst = RIBF(q, a, b, 1)})


system.time({textbookRst = RIBF(q, a, b, 2, 1e-38, 1e-13)})


range(textbookRst - Rrst)
range(boostRst - Rrst)




range(2 * (textbookRst - Rrst) / (textbookRst + Rrst), na.rm = T)


range( 2 * ( boostRst - Rrst ) / (Rrst + boostRst), na.rm = T)




  
  
  
  
  