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




#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// testMonoLinearInterpo
NumericVector testMonoLinearInterpo(NumericVector x, NumericVector y, NumericVector xtest, String monotype);
RcppExport SEXP sourceCpp_1_testMonoLinearInterpo(SEXP xSEXP, SEXP ySEXP, SEXP xtestSEXP, SEXP monotypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xtest(xtestSEXP);
    Rcpp::traits::input_parameter< String >::type monotype(monotypeSEXP);
    rcpp_result_gen = Rcpp::wrap(testMonoLinearInterpo(x, y, xtest, monotype));
    return rcpp_result_gen;
END_RCPP
}
