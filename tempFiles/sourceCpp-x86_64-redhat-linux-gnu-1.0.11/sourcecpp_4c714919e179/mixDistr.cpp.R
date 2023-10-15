`.sourceCpp_60_DLLInfo` <- dyn.load('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c714919e179/sourceCpp_61.so')

mixDistrList <- Rcpp:::sourceCppFunction(function(X, Y, rstSizes, Yweights, regridMethod = "r4", maxCore = 1000L) {}, FALSE, `.sourceCpp_60_DLLInfo`, 'sourceCpp_60_mixDistrList')
mixDistr <- Rcpp:::sourceCppFunction(function(X, Y, rstSize, yweight, regridMethod = "r4") {}, FALSE, `.sourceCpp_60_DLLInfo`, 'sourceCpp_60_mixDistr')

rm(`.sourceCpp_60_DLLInfo`)
