`.sourceCpp_32_DLLInfo` <- dyn.load('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_a12f4c0002b9/sourceCpp_69.so')

p0ensembleLDSlinear <- Rcpp:::sourceCppFunction(function(rawMdr, dr, windowSize, slidingSpeed, eps, targetMDRs, Nsample, sampleSize, seed, maxCore) {}, FALSE, `.sourceCpp_32_DLLInfo`, 'sourceCpp_32_p0ensembleLDSlinear')
makeP0test <- Rcpp:::sourceCppFunction(function(mdr, dr, windowSize, slidingSpeed, eps) {}, FALSE, `.sourceCpp_32_DLLInfo`, 'sourceCpp_32_makeP0test')
p0ensembleLDSlinearOld <- Rcpp:::sourceCppFunction(function(rawMdr, dr, windowSize, slidingSpeed, eps, targetMDRs, Nsample, sampleSizeRatio, seed, maxCore) {}, FALSE, `.sourceCpp_32_DLLInfo`, 'sourceCpp_32_p0ensembleLDSlinearOld')

rm(`.sourceCpp_32_DLLInfo`)
