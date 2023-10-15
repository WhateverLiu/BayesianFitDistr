distlist = lapply(1:10, function(x) list(val = sort(runif(10)), P = runif(10)))
for (i in 1:length(distlist)) 
  distlist[[i]]$P = distlist[[i]]$P / sum(distlist[[i]]$P)
MDR = runif(10)
rst = extractMain(distlist, MDR, normalizeMainPart = TRUE,
            removeZeroPs = FALSE)
