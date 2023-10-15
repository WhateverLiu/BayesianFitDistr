// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
// #include <R.h>
using namespace Rcpp;
#include "../h/charlieThreadPool.hpp"


// [[Rcpp::export]]
NumericVector checkPbetaThreadSafe(
    NumericVector q, NumericVector alpha, NumericVector beta, 
    int maxCore, int grainSize)
{
  
  CharlieThreadPool cp(maxCore);
  maxCore = cp.maxCore;
  double *x = &q[0], *a = &alpha[0], *b = &beta[0];
  NumericVector rstval(q.size());
  double *rst = &rstval[0];
  cp.parFor(0, q.size(), [x, a, b, rst](std::size_t i, std::size_t t)->bool
  {
    rst[i] = R::pbeta(x[i], a[i], b[i], 1, 1);
    return false;
  }, grainSize, 
  [](std::size_t t)->bool{ return false; }, 
  [](std::size_t t)->bool{ return false; });
  
  
  return rstval;
}


// [[Rcpp::export]]
NumericVector checkQbetaThreadSafe(
    NumericVector p, NumericVector alpha, NumericVector beta, 
    int maxCore, int grainSize)
{
  
  CharlieThreadPool cp(maxCore);
  maxCore = cp.maxCore;
  double *x = &p[0], *a = &alpha[0], *b = &beta[0];
  NumericVector rstval(p.size());
  double *rst = &rstval[0];
  cp.parFor(0, p.size(), [x, a, b, rst](std::size_t i, std::size_t t)->bool
  { 
    rst[i] = R::qbeta(x[i], a[i], b[i], 1, 0);
    return false;
  }, grainSize,  
  [](std::size_t t)->bool{ return false; }, 
  [](std::size_t t)->bool{ return false; });
  
  
  return rstval;
} 














