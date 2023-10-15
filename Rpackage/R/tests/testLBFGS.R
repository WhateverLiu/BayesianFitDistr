
debugRcpp.sh


Rcpp::sourceCpp("cpp/tests/testLBFGSandPlainGD.cpp")


init = runif(4) * 100
miniRosenbrock(init, 1000, giveGrad = F, lb = c(0.1, 0.1, 0.2, 0.3), ub = c(1.9, 0.9, 0.8, 1.1))
miniRosenbrock(init, 1000, giveGrad = F)

# miniRosenbrockPlainGD(init, 1000, giveGrad = T, lr = 1e-4)

