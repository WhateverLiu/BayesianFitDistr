`.sourceCpp_23_DLLInfo` <- dyn.load('/finance_develop/Charlie/BayesianFitDistr/Rpackage/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_1641452e47543/sourceCpp_24.so')

mixDistrList <- Rcpp:::sourceCppFunction(function(X, Y, rstSizes, Yweights, regridMethod = "r4", maxCore = 1000L) {}, FALSE, `.sourceCpp_23_DLLInfo`, 'sourceCpp_23_mixDistrList')
mixDistr <- Rcpp:::sourceCppFunction(function(X, Y, rstSize, yweight, regridMethod = "r4") {}, FALSE, `.sourceCpp_23_DLLInfo`, 'sourceCpp_23_mixDistr')

rm(`.sourceCpp_23_DLLInfo`)
