#include <Rcpp.h>
using namespace Rcpp;
#define vec std::vector
#include "../h/distance.hpp"
#include "../h/trbCharlie.hpp"



// [[Rcpp::export]]
double testNegllhDiscrete(DataFrame X, NumericVector abc_lm1, 
                          double eps = 1e-8, int maxit = 100)
{
  TrB trb;
  trb.reset<true> (
      abc_lm1[0], abc_lm1[1], abc_lm1[2], abc_lm1[3], eps, maxit);
  NegllhDiscrete D;
  NumericVector x = X[0], p = X[1];
  auto cdf = [&trb](double q)->double { return trb.cdf(q); };
  return D(&x[0], &p[0], x.size(), cdf);
}



#undef vec