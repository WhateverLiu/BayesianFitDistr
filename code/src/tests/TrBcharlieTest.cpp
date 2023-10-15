// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <R.h>
using namespace Rcpp;
#define vec std::vector
#include "../h/trbCharlie.hpp"



// [[Rcpp::export]]
List testTrBcharlie(NumericVector x, NumericVector abcd) 
{
  NumericVector logpdf(x.size()), cdf(x.size());
  TrB trb; 
  trb.reset(abcd[0], abcd[1], abcd[2], abcd[3]);
  for (int i = 0, iend = x.size(); i < iend; ++i)
  {
    logpdf[i] = trb.logPdf(x[i]); 
    cdf[i] = trb.cdf(x[i]);
  }
  auto l1m = trb.computeL1m(abcd[3]);
  return List::create(
    Named("logpdf") = logpdf,
    Named("cdf") = cdf,
    Named("l1m") = l1m.first
  );
}


// [[Rcpp::export]]
double testCharlieTrbReset(double a, double b, double c, double lm1, 
                           double eps, int maxit, bool useNewton)
{
  TrB trb;
  if (useNewton) return trb.reset<true>(a, b, c, lm1, eps, maxit);
  return trb.reset<false>(a, b, c, lm1, eps, maxit);
}



// double xfirst, double xlast, int size,
// double a, double b, double c, double lm1, 
// double eps, int maxit,
// double lastProbLowerBound,
// double *pmf
// Tail probability will be loaded to the last point.
// [[Rcpp::export]]
List discretize(
    NumericVector abc_lm1, 
    double lastPoint, int size, 
    double eps, int maxit, 
    double lastProbLowerBound,
    bool useNewton = false)
{ 
  size = std::max(3, size);
  TrB trb;
  double firstPoint = lastPoint / size;
  NumericVector val(size), p(size);
  if (useNewton)
  { 
    trb.discretize<true> (
        firstPoint, lastPoint, size, 
        abc_lm1[0], abc_lm1[1], abc_lm1[2], abc_lm1[3],
        eps, maxit, 
        lastProbLowerBound, &val[0], &p[0]);
  } 
  else
  { 
    trb.discretize<false> (
        firstPoint, lastPoint, size, 
        abc_lm1[0], abc_lm1[1], abc_lm1[2], abc_lm1[3],
        eps, maxit, 
        lastProbLowerBound, &val[0], &p[0]);
  } 
  // NumericVector val(size);
  double delta = firstPoint;
  for (int i = 0; i < size; ++i) val[i] = firstPoint + i * delta;
  
  
  NumericVector abcd(abc_lm1.begin(), abc_lm1.end());
  abcd[3] = trb.d;
  return List::create(
    Named("abcd") = abcd, 
    Named("distr") = DataFrame::create(Named("val") = val, Named("P") = p)
  ); 
}
















#undef vec







