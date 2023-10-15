

cd /finance_develop/Charlie/BayesianFitDistr

debugRcpp.sh


Rcpp::sourceCpp('Rpackage/src/VecPool.cpp', cacheDir = '../tempFiles', 
                verbose = 1)
testVecPool()









