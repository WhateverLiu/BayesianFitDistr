library(NGFMfitDistr)
data(tmpCondEmpDistrs, package = "NGFMfitDistr")
data(tmpLm1, package = "NGFMfitDistr")


tmpOptRst = LBFGSBtrbFitList(
  abc = matrix(c(4, 5, 6)), 
  lm1 = tmpLm1,
  empDistrList = tmpCondEmpDistrs,
  abcLB = c(1.01, 0.1, 0.1), abcUB = c(30, 30, 30), 
  scaleEps = 1e-8, scaleMaxit = 100, 
  distanceFun = "likelihood",
  max_iterations = 100, maxCore = 1000, 
  RIBlib = "Numerical Recipes", 
  sequentialUpdate = 5000, 
  hgrad = 0, centralDiff = TRUE, m = 6, epsilon = 1e-5, 
  epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
  max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
  ftol = 1e-4, wolfe = 0.9)


if (max(abs(tmpOptRst$fval)) > 1e200)
  warning("Not all have converged.")


plot(tmpOptRst$fval, xlab = "MDR index", type = "l", ylab = "Obj value", 
     las = 1, bty = "L")