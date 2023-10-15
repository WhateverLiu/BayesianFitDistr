// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../h/monoLinearInterpo.hpp"


// [[Rcpp::export]]
NumericVector testMonoLinearInterpo(NumericVector x, NumericVector y, 
                                    NumericVector xtest, 
                                    String monotype = "increasing") 
{
  NumericVector rst(xtest.size());
  MonoLinearInterpo<int, double, double> mli;
  if      (monotype == "increasing") 
    mli.reset<0> (&*x.begin(), &*x.end(), &*y.begin());
  else if (monotype == "nondecreasing")
    mli.reset<1> (&*x.begin(), &*x.end(), &*y.begin());
  else if (monotype == "decreasing")
    mli.reset<2> (&*x.begin(), &*x.end(), &*y.begin());
  else 
    mli.reset<3> (&*x.begin(), &*x.end(), &*y.begin());
  
  
  mli(&*xtest.begin(), &*xtest.end(), &*rst.begin());
  return rst;
}


