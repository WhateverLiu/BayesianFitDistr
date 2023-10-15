// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
void assignP0(NumericVector p, double p0) 
{
  double main = std::accumulate(p.begin() + 1, p.end(), 0.0);
  double r = (1.0 - p0) / main;
  p[0] = p0;
  for (auto it = p.begin() + 1; it < p.end(); ++it) *it *= r;
}


// // [[Rcpp::export]]
// void inplaceNormalize(NumericVector p)
// {
//   double r = 1.0 / std::accumulate(p.begin(), p.end(), 0.0);
//   for (auto it = p.begin(); it != p.end(); ++it) *it *= r;
// }


#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// assignP0
void assignP0(NumericVector p, double p0);
RcppExport SEXP sourceCpp_38_assignP0(SEXP pSEXP, SEXP p0SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type p0(p0SEXP);
    assignP0(p, p0);
    return R_NilValue;
END_RCPP
}
