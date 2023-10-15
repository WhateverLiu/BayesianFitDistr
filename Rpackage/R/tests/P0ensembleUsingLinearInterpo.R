


cd /finance_develop/Charlie/BayesianFitDistr
debugRcpp.sh


Rcpp::sourceCpp("Rpackage/src/P0ensembleUsingLinearInterpo.cpp", 
                cacheDir = '../tempFiles', verbose = 1)


set.seed(123)
N = 100000
nzero = round(N * 0.3)
windowSize = 137
slidingSpeed = 37
Nsample = 100
sampleSize = round(N * 0.632)
tmpMdr = round(sort(runif(N)), 3)
tmpCdr = runif(N)
tmpCdr[sample(N, nzero)] = 0
targetMDRs = c(1:10000, seq(10000L + 10L, 100000L-10L, by = 10L))/ 1e5
tmp = p0ensembleLDSlinear(
  tmpMdr, dr = tmpCdr, windowSize = windowSize, slidingSpeed = slidingSpeed,
  targetMDRs = targetMDRs, Nsample = Nsample, eps = 1e-10,
  sampleSize = round(0.632 * N), seed = 42, 
  maxCore = parallel::detectCores()
  )






plot(tmpMdr, tmpP0, xlim = c(0, 1), ylim = c(0, 1))
for (u in tmp$p0models)
  lines(x = targetMDRs, y = u)
lines(x = targetMDRs, y= tmp$P0, type = "l", col = "red")




set.seed(123)
N = 100000
tmpMdr = round(sort(runif(N)), 2)
tmpP0 = round(runif(N), 2)
targetMDRs = c(1:10000, seq(10000L + 10L, 100000L-10L, by = 10L))/ 1e5
tmp = p0ensembleLDSlinear(
  mdr = tmpMdr, P0 = tmpP0, targetMDRs = targetMDRs,
  Nsample = 100, sampleSize = round(N * 0.632), seed =42, maxCore = 1000)
plot(tmpMdr, tmpP0, xlim = c(0, 1), ylim = c(0, 1))
for (u in tmp$p0models)
  lines(x = targetMDRs, y = u)
lines(x = targetMDRs, y= tmp$P0, type = "l", col = "red")











