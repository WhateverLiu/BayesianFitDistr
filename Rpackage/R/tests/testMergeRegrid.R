

Rcpp::sourceCpp("cpp/tests/testMergeRegrid.cpp")


i = 0L
while (T) 
{
  i = i + 1L
  x = round(runif(100), 2)
  pmf = data.frame(val = sort(runif(100)), P = runif(100))
  pmf$P = pmf$P / sum(pmf$P)
  r = runif(1)
  if (r < 1/3) w = 0
  else if (r < 2/3) w = 1
  else w = runif(1)
  
  
  regridMethod = sample(c("lr", "rg4", "lmm"), 1)
  rstSize = sample(1000, 1) + 1L
  tmp = testMergeRegrid(x, pmf, w = w, rstSize = rstSize,
                        regridMethod = regridMethod)
  rgMean = mean(x) * (1 - w) + keyALGs::Mean(pmf) * w
  rgVar = mean(x ^ 2) * (1-w) + keyALGs::Mean(
    data.frame(val = pmf$val ^ 2, P = pmf$P)) * w  - rgMean ^ 2
  
  
  if (regridMethod == "lr" || regridMethod == "lmm" || 
      (regridMethod == "rg4" && rstSize < 10L))
  {
    if (abs(keyALGs::Mean(tmp) - rgMean) > 1e-10) break
  }
  else
  {
    err = abs(keyALGs::Mean(tmp) - rgMean) + abs(keyALGs::Var(tmp) - rgVar)
    if (err > 1e-10) break
  }
}




# debugRcpp.sh
source("R/sourceCppCharlie.R")
Rcpp::sourceCpp("cpp/tests/testMakeEmpDistr.cpp")
N = 100000L
Nwin = 10000L
x = round(runif(N), 4)
winSize = 1000L
tmp = sample(N - winSize, Nwin)
wins = rbind(tmp, tmp + winSize - 1L)
# wins = matrix(sample(N, Nwin * 2L), nrow = 2)
# wins = apply(wins, 2, function(x) sort(x))


# wins = wins[, order(wins[2, ] - wins[1, ], decreasing = T)]
system.time({
  myrst = makeEmpDistrList(x, wins, 50L, maxCore = 1000, regridMethod = "r4")
})
# myrstMeans= unlist(lapply(myrst, function(x) keyALGs::Mean(x)))
# myrstVars  =unlist(lapply(myrst, function(x) keyALGs::Var(x)))


system.time({trueRst = lapply(as.data.frame(wins), function(u)
{
  keyALGs::emp( x[u[1]:u[2]], N = 50L, rgMethod = "r4")
})})
range(unlist(trueRst) - unlist(myrst))




trueMeans = mapply(function(u)
{
  mean(x[u[1]:u[2]])
}, as.data.frame(wins))
trueVars = mapply(function(u)
{
  n = u[2] - u[1] + 1L
  var(x[u[1]:u[2]]) * ((n - 1) / n)
}, as.data.frame(wins))


range(myrstMeans - trueMeans)
range(myrstVars - trueVars)








source("R/sourceCppCharlie.R")
Rcpp::sourceCpp("cpp/tests/testMakeEmpDistr.cpp")
N = 10000L
Nwin = 1000L
x = round(runif(N), 4)
wins = matrix(sample(N, Nwin * 2L), nrow = 2)
wins = apply(wins, 2, function(x) sort(x))
system.time({
  distlist1 = makeEmpDistrList(x, wins, 50L, maxCore = 1000, regridMethod = "r4")
})


N = 10000L
Nwin = 1000L
x = round(runif(N), 4)
wins = matrix(sample(N, Nwin * 2L), nrow = 2)
wins = apply(wins, 2, function(x) sort(x))
system.time({
  distlist2 = makeEmpDistrList(x, wins, 50L, maxCore = 1, regridMethod = "r4")
})


Yweights = 0.3
rstSize= 20
myrst = mixDistrList(distlist1, distlist2, rstSize, Yweights = Yweights, regridMethod = "r4")
range(unlist(lapply(tmp, function(x) sum(x$P))))
# rstMeans = unlist(lapply(tmp, function(x) keyALGs::Mean(x)))
# means1 = unlist(lapply(distlist1, function(x) keyALGs::Mean(x)))
# means2 = unlist(lapply(distlist2, function(x) keyALGs::Mean(x)))
# range(rstMeans - (means1 * (1 - Yweights) + means2 * Yweights))


trueRst= mapply(function(X, Y)
{
  X[[2]] = (1 - Yweights) * X[[2]]
  Y[[2]] = Yweights * Y[[2]]
  dat = data.table::data.table(val = c(X[[1]], Y[[1]]), P = c(X[[2]], Y[[2]]))
  dat = data.table::setDF(dat[, .(P = sum(P)), by = val])
  dat = dat[order(dat$val), , drop = F]
  ng = seq(dat$val[1], dat$val[nrow(dat)], len = rstSize)
  keyALGs::rglr4split(dat, ng)
}, distlist1, distlist2, SIMPLIFY = F)


range(unlist(trueRst) - unlist(myrst))






tmp = t(data.table::setDF(data.table::fread(
  "tests/WindDamageDistribution_Com_CovA_64_20190516.csv"))[-1])


system.time({tmplist = apply(tmp, 2, function(x)
{
  list(val = seq(0, x[2], len = 64), P = x[3:length(x)])
})})


system.time({dflist = apply(tmp, 2, function(x)
{
  data.frame(val = seq(0, x[2], len = 64), P = x[3:length(x)])
})})



















