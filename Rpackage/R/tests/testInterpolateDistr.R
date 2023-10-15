

N = 10000L
k = 10000L
Nwanted = 10000L
sizeWanted = sample(100L, 1L) + 1L


X = lapply(1:N, function(i)
{
  s = sample(k, 1)
  d = list(val = sort(runif(s)), P = runif(s))
  d$P = d$P / sum(d$P); d
})
mdrs = unlist(lapply(X, function(x) keyALGs::Mean(x)))
odr = order(mdrs)
mdrs = mdrs[odr]
X = X[odr]
X = c(list(list(val = 0, P = 1)), X, list(list(val = 1, P = 1)))
mdrs = c(0, mdrs, 1)
MDRwanted = runif(Nwanted)


Rcpp::sourceCpp("cpp/tests/testInterpolateDistr.cpp")


system.time({
rst = findEmpDistrGivenMDR(X, MDR = mdrs, MDRwanted = MDRwanted, 
                     sizeWanted = sizeWanted, maxCore = 1000,
                     regridMethod = "r4")
})


# i = 0L
trueRst = lapply(MDRwanted, function(m)
{
  # i <<- i + 1L
  h = which(mdrs >= m)[1]
  l = h - 1L
  w = (m - mdrs[l]) / (mdrs[h] - mdrs[l])
  dat = data.table::data.table(
    val = c(X[[l]]$val, X[[h]]$val),
    P = c(X[[l]]$P * (1-w), X[[h]]$P * w)
  )
  dat = data.table::setDF(dat[, .(P = sum(P)), by = val])
  # print(str(dat))
  dat = dat[order(dat$val), , drop  =F]
  ng = seq(dat$val[1], dat$val[nrow(dat)], len = sizeWanted)
  keyALGs::rglr4split(dat, ng, returnDF = F)
})


range(unlist(trueRst) - unlist(rst))



































