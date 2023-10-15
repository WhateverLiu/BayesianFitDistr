`.sourceCpp_21_DLLInfo` <- dyn.load('/finance_develop/Charlie/BayesianFitDistr/Rpackage/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_164145fcccc33/sourceCpp_22.so')

makeEmpDistr <- Rcpp:::sourceCppFunction(function(x, pmf, w, rstSize, regridMethod = "r4", fixedMin = 1e300, fixedMax = 1e300, biasCorrectionMultiplier = 1.0) {}, FALSE, `.sourceCpp_21_DLLInfo`, 'sourceCpp_21_makeEmpDistr')
makeEmpDistrList <- Rcpp:::sourceCppFunction(function(X, windows, rstSizes, regridMethod = "r4", maxCore = 1000L, fixedMin = 1e300, fixedMax = 1e300, biasCorrectionMultiplier = 1.0) {}, FALSE, `.sourceCpp_21_DLLInfo`, 'sourceCpp_21_makeEmpDistrList')

rm(`.sourceCpp_21_DLLInfo`)
