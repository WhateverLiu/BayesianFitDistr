// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;


//' Extract main part.
//' 
//' Remove \eqn{P_0}s from the PMFs in a list.
//' 
//' @param distlist  A list of PMFs. Each PMF is a list of two numeric vectors.
//' The first is the support. The second is the probabilities.
//' 
//' @param MDR  A numeric vector of the prescribed MDRs. Should have the same
//' size as \code{distlist}.
//' 
//' @param normalizeMainPart  A boolean value. Should we normalize the 
//' probabilities after removing \eqn{P_0}? Default to \code{TRUE}.
//' 
//' @param removeZeroPs  A boolean value. Some probabilities can be zeros in
//' \code{distlist}. Should we remove them along with their points on the
//' support? Default to \code{FALSE}.
//' 
//' @return A list of two:
//' 
//' \code{$conditionalDistrs}:  A list of PMFs to be fitted by TrBs.
//' 
//' \code{$lm1}: A numeric vector of the target limited means 
//' for constraining \eqn{\mathbb{E}[\min(X, 1)]} of the TrB models.
//' 
//' @inherit computeWindows details
//' 
//' @example  inst/examples/extractMain.R
//' 
// [[Rcpp::export]]
List extractMain(List distlist, NumericVector MDR, 
                 bool normalizeMainPart = true,
                 bool removeZeroPs = false)
{
  List rst(distlist.size());
  NumericVector lm1(distlist.size());
  if (!removeZeroPs)
  {
    for (int i = 0, iend = distlist.size(); i < iend; ++i)
    {
      List d = distlist[i];
      NumericVector val = d[0], P = d[1];
      NumericVector valNew(val.begin() + 1, val.end());
      NumericVector Pnew(P.begin() + 1, P.end());
      double r = normalizeMainPart ? 1.0 / (1.0 - P[0]) : 1.0;
      for (auto it = Pnew.begin(); it < Pnew.end(); ++it) *it *= r;
      lm1[i] = MDR[i] * r;
      rst[i] = List::create(Named("val") = valNew, Named("P") = Pnew);
    }
  }
  else
  {
    for (int i = 0, iend = distlist.size(); i < iend; ++i)
    {
      List d = distlist[i];
      NumericVector val = d[0], P = d[1];
      int Nnonzero = 0;
      for (auto it = P.begin() + 1; it < P.end(); ++it) Nnonzero += *it != 0;
      NumericVector valNew(Nnonzero);
      NumericVector Pnew(Nnonzero);
      double r = normalizeMainPart ? 1.0 / (1.0 - P[0]) : 1.0;
      int h = 0;
      for (int k = 1, kend = P.size(); k < kend; ++k)
      {
        if (P[k] == 0) continue;
        valNew[h] = val[k];
        Pnew[h] = P[k] * r;
        ++h;
      }
      lm1[i] = MDR[i] * r;
      rst[i] = List::create(Named("val") = valNew, Named("P") = Pnew);
    } 
  }
  
  
  return List::create(Named("conditionalDistrs") = rst, Named("lm1") = lm1);
}





















#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// extractMain
List extractMain(List distlist, NumericVector MDR, bool normalizeMainPart, bool removeZeroPs);
RcppExport SEXP sourceCpp_46_extractMain(SEXP distlistSEXP, SEXP MDRSEXP, SEXP normalizeMainPartSEXP, SEXP removeZeroPsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type distlist(distlistSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MDR(MDRSEXP);
    Rcpp::traits::input_parameter< bool >::type normalizeMainPart(normalizeMainPartSEXP);
    Rcpp::traits::input_parameter< bool >::type removeZeroPs(removeZeroPsSEXP);
    rcpp_result_gen = Rcpp::wrap(extractMain(distlist, MDR, normalizeMainPart, removeZeroPs));
    return rcpp_result_gen;
END_RCPP
}
