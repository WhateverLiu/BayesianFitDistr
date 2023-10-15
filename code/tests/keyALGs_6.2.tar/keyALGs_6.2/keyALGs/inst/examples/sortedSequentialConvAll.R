

# Generate a list of 100 PMFs.
distlist = lapply(1:100, function(i)
{
  supportSize = sample(200, 1)
  supportRange = sort(runif(2))
  support = seq(supportRange[1], supportRange[2], len = supportSize)
  X = data.frame(val = support, P = runif(supportSize))
  X$P = X$P / sum(X$P); X
})


seqAggDistr = keyALGs::sequentialConvAll(distlist)
sortedSeqAggDistr = keyALGs::sortedSequentialConvAll(distlist)


theoreticalVar = sum(unlist(lapply(distlist, function(x) keyALGs::Var(x))))
cat("Theoretical aggregate variance = ", theoreticalVar, "\n",
    "Sequential aggregate variance = ", keyALGs::Var(seqAggDistr), "\n",
    "Sorted sequential aggregate variance = ", 
    keyALGs::Var(sortedSeqAggDistr), "\n", 
    "Sorted sequential accumulation MAY mitigate variance inflation.",
    sep = "")



















