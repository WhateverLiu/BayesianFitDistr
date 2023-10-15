// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../h/optimizer.hpp"


// [[Rcpp::export]]
IntegerMatrix testDeltas(int dim) 
{
  Deltas D(dim);
  IntegerMatrix rst(D.vecDim(), D.Nvec());
  std::copy(D.v.begin(), D.v.end(), rst.begin());
  return rst;
}


struct Rosen
{
  double operator()(double *p)
  {
    double x = p[0], y = p[1];
    double t = y - x * x;
    return (1 - x) * (1 - x) + 100 * t * t;
  }
};


struct LR
{
  void operator()(double *x, double *lr)
  {
    lr[0] = 0.0001;
    lr[1] = 0.0001;
  }
};


// [[Rcpp::export]]
NumericVector testGridSearch(NumericVector xinit, int maxit)
{
  GridSearch<2, int16_t> G(maxit);
  LR lr;
  Rosen f;
  NumericVector rst(xinit.begin(), xinit.end());
  // (double *x, Fun &f, LR &lr)
  G(&rst[0], f, lr);
  return rst;
}



















