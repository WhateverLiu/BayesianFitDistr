

source('R/sourceCppCharlie.R')
sourceCppCharlie('cpp/tests/bisecRoot.cpp', verbose = T)


bisecRootTest(runif(1), eps = 1e-6)
bisecRootTest(runif(1) * 10, eps = 1e-6)









