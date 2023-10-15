// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#define vec std::vector


// constexpr const bool useBoost = false;
#include <boost/math/special_functions/beta.hpp>
// [[Rcpp::depends(BH)]]
#include "../h/textbookRegularizedIncompleteBeta.h"
#include "../h/regularizedIncompleteBetaFun.h"
#include "../h/trbCharlie.hpp"
#include "../solve-d.hpp"


// [[Rcpp::export]]
List testLimitedMean(NumericVector abc_lm1, double limit, int maxit = 100, 
                       double eps = 1e-8) 
{
  TrB trb;
  trb.reset<true, 2> (abc_lm1[0], abc_lm1[1], abc_lm1[2],
                  abc_lm1[3], eps, maxit);
  double lev = trb.computeLimitedMean<2> (limit);
  return List::create(Named("d") = trb.d, Named("limitedMean") = lev);
}













