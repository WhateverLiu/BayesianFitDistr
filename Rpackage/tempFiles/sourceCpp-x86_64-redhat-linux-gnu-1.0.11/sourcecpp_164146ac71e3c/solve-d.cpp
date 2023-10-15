// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "h/trbCharlie.hpp"




template <int RIBlib>
NumericVector solve_d(NumericVector abc_lm1, double eps, int maxit, 
                      bool useNewton) 
{
  TrB trb;
  double *param = &abc_lm1[0];
  NumericVector rst(abc_lm1.size() / 4);
  for (int i = 0, iend = abc_lm1.size(); i < iend; i += 4)
  {
    double *p = param + i;
    if (useNewton) trb.reset<true, RIBlib> (p[0], p[1], p[2], p[3], eps, maxit);
    else trb.reset<false, RIBlib> (p[0], p[1], p[2], p[3], eps, maxit);
    rst[i >> 2] = trb.d;
  }
  return rst;
}




//' Solve the scale parameter d.
//' 
//' Solve TrB parameter d given a, b, c and the constrained limited mean at 1.
//' 
//' @param abc_lm1  A numeric matrix with four rows. Each column is a vector of
//' \code{(a, b, c)} and the target limited mean at 1.
//' 
//' @param eps  If the difference between TrB's limited mean and \code{lm1} is
//' less than \code{eps}, stop the Newton's iteration. Default to \code{1e-8}.
//' 
//' @param maxit  The maximum number of iterations for Newton's method.
//' Default to \code{100}.
//' 
//' @param useNewton  A boolean. \code{TRUE} invokes Newton's method. 
//' \code{FALSE} invokes bisection. Default to \code{TRUE}.
//' 
//' @inheritParams LBFGSBtrbFitList
//' 
//' @return A numeric vector of \code{d}s. The vector is of of size \code{ncol(abc_lm1)}. 
//' 
// [[Rcpp::export]]
NumericVector solve_d(NumericVector abc_lm1, double eps = 1e-8, int maxit = 100, 
                      bool useNewton = true, String RIBlib = "R::pbeta")
{
  if (RIBlib == "Boost") return solve_d<0> (abc_lm1, eps, maxit, useNewton);
  if (RIBlib == "R::pbeta") return solve_d<1> (abc_lm1, eps, maxit, useNewton);
  if (RIBlib == "Numerical Recipes") return solve_d<2> (abc_lm1, eps, maxit, useNewton);
  stop("Method of computing Regularized Incomplete Beta function is not implemented.");
  return NumericVector(0);
}
  






#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// solve_d
NumericVector solve_d(NumericVector abc_lm1, double eps, int maxit, bool useNewton, String RIBlib);
RcppExport SEXP sourceCpp_27_solve_d(SEXP abc_lm1SEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP useNewtonSEXP, SEXP RIBlibSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type abc_lm1(abc_lm1SEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< bool >::type useNewton(useNewtonSEXP);
    Rcpp::traits::input_parameter< String >::type RIBlib(RIBlibSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_d(abc_lm1, eps, maxit, useNewton, RIBlib));
    return rcpp_result_gen;
END_RCPP
}
