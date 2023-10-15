// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#define vec std::vector


// constexpr const bool useBoost = false;
#include <boost/math/special_functions/beta.hpp>
// [[Rcpp::depends(BH)]]


constexpr const bool useBoost = true;
#include "../h/regularizedIncompleteBetaFun.h"


// Regularized incomplete beta function.
// [[Rcpp::export]]
NumericVector RIBF(NumericVector q, NumericVector a, NumericVector b, 
                   bool useboost = false)
{
  if (a.size() != b.size()) stop("a and b have different sizes.");
  int asize = a.size();
  NumericVector rst(q.size());
  for (int i = 0, iend = q.size(); i < iend; ++i)
  {
    int k = i % asize;
    rst[i] = ribetaf(q[i], a[k], b[k], 1, 0, useboost);
  }
  return rst;
}
















































