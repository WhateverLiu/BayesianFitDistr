# Generate a list of 100 PMFs.
distlist = lapply(1:100, function(i)
{
  supportSize = sample(200, 1)
  supportRange = sort(runif(2))
  support = seq(supportRange[1], supportRange[2], len = supportSize)
  X = data.frame(val = support, P = runif(supportSize))
  X$P = X$P / sum(X$P); X
})
aggDistr = keyALGs::sequentialConvAll(distlist)
plot(aggDistr, type = "h")
