

source("R/rfuns.R")
source("R/sourceCppCharlie.R")
sourceCppCharlie("cpp/amalgamation.cpp", Oflag = '-O3', verbose = T)


abc = 10 ^ runif(3, 0, 1)
abc_lm1 = c(abc, runif(1))
d = solve_d(abc_lm1, eps = 1e-8, 100, useNewton = T) 
abcd = c(abc, d)


rs = actuar_rtrbeta(1000, abcd)
rs = pmin(1, rs)
distr = keyALGs::emp(rs, 20, "r4")


init = sort(10 ^ runif(3, 0, 1.5))
init[1:2] = sample(init[1:2])
init = c(init, abc_lm1[4])


rst = LBFGSBtrbFit(
  init, empDistr = distr, 
  abcLB = c(1.001, 0.001, 0.001), 
  abcUB = c(30, 30, 30),
  scaleEps = 1e-8, scaleMaxit = 100, distanceFun = "likelihood",
  max_iterations = 100, hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5,
  epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10,
  max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, ftol = 1e-4,
  wolfe = 0.9); rst; abc_lm1; init


optabcd = c(rst$param[1:3], solve_d(rst$param, 1e-8, 100))
optDistr = actuar_discretize(optabcd, distr$val)
optDistrH = wrapAsHist(optDistr)
ylim = range(c(0, max(distr$P), max(optDistr$P)))
plot(optDistrH, col = "gray", freq = F, border = NA, ylim = ylim, 
     bty = "L", las = 1)
lines(distr, type = "h", lwd = 3)












source("R/rfuns.R")
source("R/sourceCppCharlie.R")
Rcpp::sourceCpp("cpp/amalgamation.cpp")
# sourceCppCharlie("cpp/amalgamation.cpp", Oflag = '-O3', verbose = T)


set.seed(123)


# N = 1
# a = runif(N, min = log(1.001), max = log(30))
# b = runif(N, min = log(0.01), max = log(30))
# c = runif(N, min = log(1), max = log(30))
# lm1 = runif(N, min = log(1e-10), max = log(1))
a = log(2); b = log(3); c = log(4); lm1 = log(0.5)
s = 0.01
N = 1000L
lb = c(log(1.001), log(0.01), log(1), log(1e-10))
ub = c(log(30), log(30), log(30), log(1))
params = list(c(a, b, c, lm1))
for ( i in 2:N)
{
  tmp = params[[i - 1]] + rnorm(4, 0, s)
  params[[i]] = pmax(lb, pmin(tmp, ub))
}
abc = exp(matrix(unlist(params), nrow = 4))
lm1 = abc[4, , drop = T]
abc = abc[1:3, , drop = F]


abcLB = exp(lb)
abcUB = exp(ub)


empDistrs = mapply(function(abc, lm1)
{
  d = solve_d(c(abc, lm1), 1e-8, 100)
  x = actuar_rtrbeta(1000, abcd = c(abc, d))
  x = pmin(x, 1)
  keyALGs::emp(x, N = 20, rgMethod = "r4")
}, as.data.frame(abc), lm1, SIMPLIFY = F)


abc = exp(seq(log(1.001), log(30), len = 5)[2:4])
abc = t(expand.grid(a = abc, b = abc, c = abc))
lm1 = unlist(lapply(empDistrs, function(x) keyALGs::Mean(x)))
names(lm1) = NULL


save(abc, lm1, empDistrs, abcLB, abcUB, file = "tests/tmpData.Rdata")


# ==============================================================================
# Debug session.
# ==============================================================================
source("R/rfuns.R")
source("R/sourceCppCharlie.R")
Rcpp::sourceCpp("cpp/amalgamation.cpp")
load("tests/tmpData.Rdata")
# ==============================================================================
system.time({
optRst = LBFGSBtrbFitList(
  abc = abc, lm1 = lm1, empDistrList = empDistrs, abcLB = abcLB, abcUB = abcUB, 
  scaleEps = 1e-8, scaleMaxit = 100, distanceFun = "likelihood", 
  max_iterations = 100, maxCore = 1000, RIBlib = "Numerical Recipes", 
  sequentialUpdate = F, hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
  epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
  max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
  ftol = 1e-4, wolfe = 0.9)
})


# tmp= colSums(-optRst$param[1:3, , drop = F] * log(abctrue))
# odr = order(tmp, decreasing = T)
# optRst$param[1:3, odr[1:5], drop = F]
# abctrue[, odr[1:5], drop = F]


pdf("tests/test.pdf", width = 10, height = 10 * (9/16))
par(mar = c(0, 0, 0, 0), mfrow = c(4, 5))
i = 0L
tmp = mapply(function(x, distr)
{
  i <<- i + 1L
  optabcd = c(x[1:3], solve_d(x, 1e-8, 100))
  optDistr = actuar_discretize(optabcd, distr$val)
  distrH = wrapAsHist(distr)
  ylim = range(c(0, max(distr$P), max(optDistr$P)))
  plot(distrH, col = "gray", freq = F, border = NA, ylim = ylim, 
       bty = "L", las = 1, xaxt  ="n", yaxt = "n", xlab = "", ylab = "",
       main = "")
  wd = diff(distrH$mids[1:2])
  wd = wd / 8
  rect(xleft = optDistr$val - wd, ybottom = 0, 
       xright = optDistr$val + wd, ytop = optDistr$P,
       col = "black", border = NA)
  legend("top", legend = i, bty = "n", cex = 2)
}, as.data.frame(optRst$param), empDistrs, SIMPLIFY = F)
dev.off()




rst = LBFGSBtrbFit(
  # c(abc[, 1], lm1[1]),
  c(abc[, 1], lm1[1]),
  # c(2,3, 4,0.5),
  empDistr = empDistrs[[1]], 
  abcLB = c(1.001, 0.001, 0.1), 
  abcUB = c(30, 30, 30),
  scaleEps = 1e-8, scaleMaxit = 100, distanceFun = "likelihood",
  max_iterations = 100, hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5,
  epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10,
  max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, ftol = 1e-4,
  wolfe = 0.9); rst; abc_lm1; init




# Why do we prefer brute-force:
# The same empirical PMF could be fitted by largely different sets of parameters.
# Observations: likelihood criterion is much superior to Komogorov distance.


































