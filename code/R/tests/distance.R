source("R/rfuns.R")
Rcpp::sourceCpp('cpp/amalgamation.cpp')


l = 1
u = 20
abc = sort(runif(3, l, u))
abc[1:3] = sort(abc[1:3])
abc[1:2] = sample(abc[1:2])
lm1 = runif(1)
abc_lm1 = c(abc, lm1)
d = solve_d(abc_lm1, 1e-8, 100, T)
abcd = c(abc, d)
actuar_levtrbeta(abcd); lm1


sp = actuar_rtrbeta(1000, abcd)
sp = pmin(sp, 1)
drs = keyALGs::emp(sp, 15, "r4")
drsH = wrapAsHist(drs)
plot(drsH, freq = F, las = 1)


Rcpp::sourceCpp("cpp/tests/distance.cpp")
myRst = testNegllhDiscrete(
  X = drs, abc_lm1 = abc_lm1, eps = 1e-8, maxit = 100)
dt = (drs$val[2] - drs$val[1]) / 2
sp = c(drs$val - dt, drs$val[nrow(drs)] + dt)
if (drs$val[nrow(drs)] - 1 > -1e-10)
{
  cf = actuar_ptrbeta(sp, abcd)
  cf[length(cf)] = 1
} else
{
  cf = actuar_ptrbeta(sp, abcd)
}
discretized = data.frame(val = drs$val, P = diff(cf))
discretized$P = discretized$P / sum(discretized$P)
actuarRst = -mean(log(discretized$P) * drs$P)
actuarRst; myRst























