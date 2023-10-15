

Rcpp::sourceCpp("cpp/tests/testLimitedMean.cpp")


abc_lm1 = runif(4)
abc_lm1[1] = runif(1, 1, 3)
x = runif(1)


rst = testLimitedMean(abc_lm1, x, maxit = 100, eps = 1e-8) 
source("R/rfuns.R")
rst$limitedMean
levtrbeta(abc_lm1, limit = x, eps = 1e-8, maxit = 100)
actuar::levtrbeta(x, abc_lm1[1] / abc_lm1[3], abc_lm1[3], 
                  abc_lm1[2] / abc_lm1[3], scale = rst$d)





