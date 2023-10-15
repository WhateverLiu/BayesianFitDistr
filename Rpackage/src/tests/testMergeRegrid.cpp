// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../h/regrid.h"
#include "../h/mergeRegrid.h"


// [[Rcpp::export]]
DataFrame testMergeRegrid(NumericVector x, DataFrame pmf, double w, 
                          int rstSize, String regridMethod = "rg4")
{
  
  MergeRegrid mr;
  NumericVector y = pmf[0], yp = pmf[1];
  NumericVector rst(rstSize), rstP(rstSize);
  if (regridMethod == "lr")
  {
    mr.operator()<0, int, double, double, double, double, double>(
        &x[0], x.size(), &y[0], &yp[0], y.size(),
        &rst[0], &rstP[0], rstSize, w);
  }
  else if (regridMethod == "lmm")
  {
    mr.operator()<1, int, double, double, double, double, double>(
        &x[0], x.size(), &y[0], &yp[0], y.size(),
        &rst[0], &rstP[0], rstSize, w);
  }
  else if (regridMethod == "rg4")
  {
    mr.operator()<2, int, double, double, double, double, double>(
        &x[0], x.size(), &y[0], &yp[0], y.size(),
        &rst[0], &rstP[0], rstSize, w);
  }
  else stop("Regrid method not implemented.");
  
  
  return DataFrame::create(Named("val") = rst, Named("P") = rstP);
}





































