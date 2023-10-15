Rdevel


source("R/rfuns.R")
Rcpp::sourceCpp('cpp/amalgamation.cpp', verbose = F)
l = 1
u = 20


abc = sort(runif(3, l, u))
abc[1:3] = sort(abc[1:3])
abc[1:2] = sample(abc[1:2])
lm1 = runif(1)
abc_lm1 = c(abc, lm1)
tentMax  = 1
myRst = discretize(abc_lm1, tentMax, 
                   lastProbLowerBound = 1e-6,
                   size = 63, eps = 1e-12, maxit = 100, 
                   useNewton = F)
abcd = myRst$abcd
myRst = myRst$distr
keyALGs::Mean(myRst) - lm1


myRst = discretize(abc_lm1, tentMax, 
                   lastProbLowerBound = 1e-6,
                   size = 63, eps = 1e-12, maxit = 100, 
                   useNewton = T)
abcd = myRst$abcd
myRst = myRst$distr
keyALGs::Mean(myRst) - lm1
lm1




source("R/rfuns.R")
Rcpp::sourceCpp('cpp/amalgamation.cpp', verbose = F)


dat = data.table::fread('tests/Coverage_A_commercial.csv')
data.table::setDF(dat)
tmp = t(dat[c('a', 'b', 'c', 'm')])


abc_lm1 = tmp[,1]
abc_lm1[4] = abc_lm1[4] / (1 - dat$p0[1])
distr = discretize(abc_lm1 = abc_lm1, 
                   tentativeMax = 1,
                   lastProbLowerBound = 1e-3,
                   size = 63, 
                   eps = 1e-12, 
                   maxit = 100,
                   useNewton = F)


































source("R/rfuns.R")
Rcpp::sourceCpp('cpp/amalgamation.cpp', verbose = F)
dat = data.table::fread('tests/Coverage_A_commercial.csv')
data.table::setDF(dat)
tmp = t(dat[c('a', 'b', 'c', 'm')])


ind = 1:ncol(tmp)
# ind = 1:10
P0 = dat$p0[ind]
system.time({distrTable = makeDistrTable(
  tmp[, ind, drop = F],
  tentativeMax = 1,
  lastProbLowerBound = 1e-5,
  P0 = P0,
  size = 64, 
  eps = 1e-11, 
  maxit = 100,  
  maxCore = 1000,
  useNewton = F, RIBlib = "Numerical Recipes")})



abc_lm1_ = c(tmp[1:3, ind, drop = F], tmp[4, ind, drop = F]/(1 - P0))


p  = 1e-2
qtrbeta(1-p, abc_lm1_ )


distr = discretize(abc_lm1 = abc_lm1_, 
                   tentativeMax = 1,
                   lastProbLowerBound = p,
                   size = 63, 
                   eps = 1e-11, 
                   maxit = 100,
                   useNewton = F, RIBlib = "R::pbeta")
distr = list(val = c(0, distr$distr$val), P = c(P0, distr$distr$P))
distr$P[-1] = distr$P[-1] * (1 -P0)

  
  
  




