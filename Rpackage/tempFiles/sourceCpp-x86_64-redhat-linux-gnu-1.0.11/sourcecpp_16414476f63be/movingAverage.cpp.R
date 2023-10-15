`.sourceCpp_25_DLLInfo` <- dyn.load('/finance_develop/Charlie/BayesianFitDistr/Rpackage/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_16414476f63be/sourceCpp_26.so')

movingAverage <- Rcpp:::sourceCppFunction(function(x, windowSize, speed, returnWindows = TRUE) {}, FALSE, `.sourceCpp_25_DLLInfo`, 'sourceCpp_25_movingAverage')
movingAverageSmoothing <- Rcpp:::sourceCppFunction(function(x, windowSize, iterations) {}, FALSE, `.sourceCpp_25_DLLInfo`, 'sourceCpp_25_movingAverageSmoothing')
windowVariances <- Rcpp:::sourceCppFunction(function(x, windowSizePercentageIncrement = 0.01, maxCore = 1000L) {}, FALSE, `.sourceCpp_25_DLLInfo`, 'sourceCpp_25_windowVariances')

rm(`.sourceCpp_25_DLLInfo`)
