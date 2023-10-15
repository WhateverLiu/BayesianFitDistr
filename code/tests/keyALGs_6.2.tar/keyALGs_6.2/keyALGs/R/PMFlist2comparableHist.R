

PMFlist2comparableHist = function(distlist, PminmaxRatio = 1e-3)
{
  for(i in 1:length(distlist)) names(distlist[[i]]) = c("val", "P")
  if(min(unlist(lapply(distlist, function(x) nrow(x)))) <= 1)
    stop("distlist contains degenerated distributions.")


  xlim = range(unlist(lapply(distlist, function(x)
  {
    range(x[x$P >= max(x$P) * PminmaxRatio, 1])
  })))


  d = unlist(lapply(distlist, function(x) diff(x$val[1:2])))
  d = d / mean(d)


  distlist = mapply(function(x, y)
  {
    x = x[x$val >= xlim[1] & x$val <= xlim[2], ]
    x$P = x$P / y; x
  }, distlist, d, SIMPLIFY = F)



  lapply(distlist, function(X)
  {
    d = diff(X$val[1:2])
    breaks = c(X$val - d / 2, X$val[nrow(X)] + d / 2)
    h = hist(X$val, breaks = breaks, plot = F)
    h$counts = X$P; h
  })

}


