

Rcpp::sourceCpp("src/tests/testQAdiscretize.cpp")
source("R/rfuns.R")


param = data.table::setDF(data.table::fread("tests/Coverage_A_commercial.csv"))
param = t(param[c("m", "max", "p0", "a", "b", "c", "d")])
rst = FTDFQAcpp(
  param, supportSize = 42, regridMethod = "lmm",
  RIBlib = "Numerical Recipes",
  outputProbsInRows = F, fineDiscretizationSize = 2000, maxCore = 1)
checkDistTable(rst)


fs = c("Coverage_A_residential.csv", "Coverage_B_commercial.csv", 
       "Coverage_B_residential.csv", "Coverage_C.csv", 
       "Coverage_D_commercial_1095.csv", "Coverage_A_commercial.csv")


for (i in 1:length(fs))
{
  
  
  pth = paste0("tests/", fs[i])
  param = data.table::setDF(data.table::fread(pth))
  
  
  if (i != 5)
  {
    param$condm = param$m / (1 - param$p0)
    ds = solve_d(t(param[c("a", "b", "c", "condm")]), eps = 1e-8, maxit = 100)
    param$dnew = ds
  } else
  {
    param$dnew = param$d
  }
  
  
  if (i != 5)
  {
    maxes = pmin(1, actuar_qtrbeta(
      1 - 1e-6 / (1 - param$p0), t(param[c("a", "b", "c", "dnew")])))
    if (min(diff(maxes)) < 0)
    {
      subInd = longestNonDecreasingSubseq(maxes)
      sfun = splinefun(x = param$m[subInd], y = maxes[subInd])
      maxes = sfun(param$m)
    }
    param$maxnew = maxes
  } else
  {
    param$maxnew = param$max
  }
  inputParam = t(param[c("m", "maxnew", "p0", "a", "b", "c", "dnew")])
  
  
  rst = FTDFQAcpp(
    inputParam, supportSize = 64, regridMethod = "r4", 
    RIBlib = "Numerical Recipes",
    outputProbsInRows = F, fineDiscretizationSize = 2000, maxCore = 1000)
  
  
  cat(checkDistTable(rst), "\n")
  
  
}


exind = -c(1,ncol(rst))
plot(x = rst[1,exind], y = rst[2,exind], type = "l", col = "blue")
lines(x = rst[1,exind], y = param$max, col = "red")


plot(x = rst[1,exind], y = rst[3,exind], type = "l", col = "blue")
lines(x = rst[1,exind], y = param$p0, col = "red")


plot(x = rst[1,exind], y = rst[66,exind], type = "l", col = "blue")
lines(x = rst[1,exind], y = param$p1, col = "red")


plot(x = rst[1,exind], y = tmp[3,-1], type = "l", col = "blue")
lines(x = rst[1,exind], y = param$p0, col = "red")








