`.sourceCpp_62_DLLInfo` <- dyn.load('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c71fb4c52d/sourceCpp_63.so')

movingAverage <- Rcpp:::sourceCppFunction(function(x, windowSize, speed, returnWindows = TRUE) {}, FALSE, `.sourceCpp_62_DLLInfo`, 'sourceCpp_62_movingAverage')
movingAverageSmoothing <- Rcpp:::sourceCppFunction(function(x, windowSize, iterations) {}, FALSE, `.sourceCpp_62_DLLInfo`, 'sourceCpp_62_movingAverageSmoothing')
windowVariances <- Rcpp:::sourceCppFunction(function(x, windowSizePercentageIncrement = 0.01, maxCore = 1000L) {}, FALSE, `.sourceCpp_62_DLLInfo`, 'sourceCpp_62_windowVariances')

rm(`.sourceCpp_62_DLLInfo`)
