`.sourceCpp_52_DLLInfo` <- dyn.load('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c712ce13731/sourceCpp_53.so')

igscore <- Rcpp:::sourceCppFunction(function(X, x, redistributePwithin = TRUE) {}, FALSE, `.sourceCpp_52_DLLInfo`, 'sourceCpp_52_igscore')
igscoreMean <- Rcpp:::sourceCppFunction(function(X, truePMF, redistributePwithin = TRUE, symmetricScore = FALSE) {}, FALSE, `.sourceCpp_52_DLLInfo`, 'sourceCpp_52_igscoreMean')

rm(`.sourceCpp_52_DLLInfo`)
