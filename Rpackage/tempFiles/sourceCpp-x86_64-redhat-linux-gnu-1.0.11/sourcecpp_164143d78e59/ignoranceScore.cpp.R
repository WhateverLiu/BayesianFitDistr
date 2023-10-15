`.sourceCpp_15_DLLInfo` <- dyn.load('/finance_develop/Charlie/BayesianFitDistr/Rpackage/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_164143d78e59/sourceCpp_16.so')

igscore <- Rcpp:::sourceCppFunction(function(X, x, redistributePwithin = TRUE) {}, FALSE, `.sourceCpp_15_DLLInfo`, 'sourceCpp_15_igscore')
igscoreMean <- Rcpp:::sourceCppFunction(function(X, truePMF, redistributePwithin = TRUE, symmetricScore = FALSE) {}, FALSE, `.sourceCpp_15_DLLInfo`, 'sourceCpp_15_igscoreMean')

rm(`.sourceCpp_15_DLLInfo`)
