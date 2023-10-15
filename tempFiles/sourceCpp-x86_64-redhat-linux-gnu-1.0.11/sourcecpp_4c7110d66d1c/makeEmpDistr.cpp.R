`.sourceCpp_58_DLLInfo` <- dyn.load('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c7110d66d1c/sourceCpp_59.so')

makeEmpDistr <- Rcpp:::sourceCppFunction(function(x, pmf, w, rstSize, regridMethod = "r4", fixedMin = 1e300, fixedMax = 1e300, biasCorrectionMultiplier = 1.0) {}, FALSE, `.sourceCpp_58_DLLInfo`, 'sourceCpp_58_makeEmpDistr')
makeEmpDistrList <- Rcpp:::sourceCppFunction(function(X, windows, rstSizes, regridMethod = "r4", maxCore = 1000L, fixedMin = 1e300, fixedMax = 1e300, biasCorrectionMultiplier = 1.0) {}, FALSE, `.sourceCpp_58_DLLInfo`, 'sourceCpp_58_makeEmpDistrList')

rm(`.sourceCpp_58_DLLInfo`)
