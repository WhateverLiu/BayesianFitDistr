# Generate a list of 100 PMFs.
distlist = lapply(1:100, function(i)
{
  supportSize = sample(200, 1)
  supportRange = sort(runif(2))
  support = seq(supportRange[1], supportRange[2], len = supportSize)
  X = data.frame(val = support, P = runif(supportSize))
  X$P = X$P / sum(X$P); X
})


comoAggDist = keyALGs::comoDistList(distlist, regridMethod = "r4")
allSDs = unlist(lapply(distlist, function(x) sqrt(keyALGs::Var(x))))
aggVarIfAllPerfectlyCorrelated = sum(allSDs) ^ 2
comoAggDistVar = keyALGs::Var(comoAggDist)


theoreticalMean = sum(unlist(lapply(distlist, function(x) keyALGs::Mean(x))))
cat("Theoretical aggregate mean = ", theoreticalMean, "\n",
    "Comonotonic aggregate mean = ", keyALGs::Mean(comoAggDist), "\n",
    "Variance upper bound of the aggregate PMF (if all random variables are
perfectly correlated) = ", aggVarIfAllPerfectlyCorrelated, "\n",
    "Comonotonic aggregate variance = ", comoAggDistVar, "\n",
    "Comonotonic aggregate variance may exceed the theoretical upper bound due
to variance inflation from regriding.", 
    sep = "")
