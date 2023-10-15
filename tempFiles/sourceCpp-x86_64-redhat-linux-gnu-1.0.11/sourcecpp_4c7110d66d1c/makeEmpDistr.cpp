// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "h/mergeRegrid.h"
#include "h/charlieThreadPool2.hpp"


// [[Rcpp::export]]
DataFrame makeEmpDistr(
    NumericVector x, DataFrame pmf, double w, 
    int rstSize, String regridMethod = "r4", 
    double fixedMin = 1e300, double fixedMax = 1e300,
    double biasCorrectionMultiplier = 1.0)
{
  
  MergeRegrid mr;
  NumericVector y = pmf[0], yp = pmf[1];
  NumericVector rst(rstSize), rstP(rstSize);
  if (regridMethod == "lr")
  {
    mr.operator()<0, int, double, double, double, double, double, double>(
        &x[0], x.size(), &y[0], &yp[0], y.size(),
        &rst[0], &rstP[0], rstSize, w, fixedMin, fixedMax, biasCorrectionMultiplier);
  }
  else if (regridMethod == "lmm")
  {
    mr.operator()<1, int, double, double, double, double, double, double>(
        &x[0], x.size(), &y[0], &yp[0], y.size(),
        &rst[0], &rstP[0], rstSize, w, fixedMin, fixedMax, biasCorrectionMultiplier);
  }
  else if (regridMethod == "r4")
  {
    mr.operator()<2, int, double, double, double, double, double, double>(
        &x[0], x.size(), &y[0], &yp[0], y.size(),
        &rst[0], &rstP[0], rstSize, w, fixedMin, fixedMax, biasCorrectionMultiplier);
  }
  else stop("Regrid method not implemented.");
  
  
  return List::create(Named("val") = rst, Named("P") = rstP);
}




// window is a 2-row integer matrix. The first row is the 1-based index
// of the first element in the window over X. The second row is the 1-based 
// index of the last element in the window over X.
// w is the weight on PMFs.
template <int regridMethod>
List makeEmpDistrList(NumericVector X, 
                      IntegerMatrix windows,
                      IntegerVector rstSizes, 
                      int maxCore, double fixedMin, double fixedMax,
                      NumericVector biasCorrectionMultiplier)
{
  const int Nwindow = windows.ncol();
  if (rstSizes.size() != 1 and rstSizes.size() != Nwindow) 
    stop("rstSizes.size() != windows.ncol().");
  
  
  for (int u: rstSizes) 
  { 
    if (u < 2) stop("rstSizes has element less than 2."); 
  }
  
  
  int *win = &windows[0];
  double *x = &X[0];
  CharlieThreadPool cp(std::move(maxCore));
  maxCore = cp.maxCore;
  std::vector<MergeRegrid> MR(maxCore);
  List rst(Nwindow);
  int rstSze = rstSizes.size();
  
  
  std::vector<double*> ptrs(Nwindow * 2);
  double **valPtrs = ptrs.data();
  double **Pptrs = valPtrs + Nwindow;
  int *distrSizes = &rstSizes[0];
  for (int i = 0; i < Nwindow; ++i)
  {
    NumericVector val(distrSizes[i % rstSze]);
    valPtrs[i] = &val[0];
    NumericVector P(distrSizes[i % rstSze]);
    Pptrs[i] = &P[0];
    rst[i] = List::create(Named("val") = val, Named("P") = P);
  }
  
  
  auto f = [&](
    std::size_t i, std::size_t t)->bool
  {
    MergeRegrid &mr = MR[t];
    double *xbegin = win[2 * i] - 1 + x;
    int xsize = win[2 * i + 1] - win[2 * i] + 1;
    int bmsize = biasCorrectionMultiplier.size();
    mr.operator() <regridMethod, int, double, double, double, double, double, double> (
        xbegin, xsize, nullptr, nullptr, 0, 
        valPtrs[i], Pptrs[i], distrSizes[i % rstSze], 0.0, fixedMin, fixedMax,
        biasCorrectionMultiplier[i % bmsize]
    );
    return false;
  };
  
  
  int grainSize = Nwindow / (maxCore * maxCore) + 1;
  cp.parFor(0, Nwindow, f, grainSize, 
            [](std::size_t t)->bool { return false; },
            [](std::size_t t)->bool { return false; });
  
  
  return rst;
}




// [[Rcpp::export]]
List makeEmpDistrList(NumericVector X, 
                      IntegerMatrix windows,
                      IntegerVector rstSizes, 
                      String regridMethod = "r4",
                      int maxCore = 1000, 
                      double fixedMin = 1e300,
                      double fixedMax = 1e300,
                      NumericVector biasCorrectionMultiplier = 1.0)
{
  
  if (regridMethod == "lr")
    return makeEmpDistrList<0>(X, windows, rstSizes, maxCore, fixedMin, fixedMax,
                               biasCorrectionMultiplier);
  
  
  if (regridMethod == "lmm")
    return makeEmpDistrList<1>(X, windows, rstSizes, maxCore, fixedMin, fixedMax,
                               biasCorrectionMultiplier);
  
  
  if (regridMethod == "r4")
    return makeEmpDistrList<2>(X, windows, rstSizes, maxCore, fixedMin, fixedMax,
                               biasCorrectionMultiplier);
  
  
  stop("Regrid method not implemented");
  return List::create();
}
  











#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// makeEmpDistr
DataFrame makeEmpDistr(NumericVector x, DataFrame pmf, double w, int rstSize, String regridMethod, double fixedMin, double fixedMax, double biasCorrectionMultiplier);
RcppExport SEXP sourceCpp_58_makeEmpDistr(SEXP xSEXP, SEXP pmfSEXP, SEXP wSEXP, SEXP rstSizeSEXP, SEXP regridMethodSEXP, SEXP fixedMinSEXP, SEXP fixedMaxSEXP, SEXP biasCorrectionMultiplierSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type pmf(pmfSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type rstSize(rstSizeSEXP);
    Rcpp::traits::input_parameter< String >::type regridMethod(regridMethodSEXP);
    Rcpp::traits::input_parameter< double >::type fixedMin(fixedMinSEXP);
    Rcpp::traits::input_parameter< double >::type fixedMax(fixedMaxSEXP);
    Rcpp::traits::input_parameter< double >::type biasCorrectionMultiplier(biasCorrectionMultiplierSEXP);
    rcpp_result_gen = Rcpp::wrap(makeEmpDistr(x, pmf, w, rstSize, regridMethod, fixedMin, fixedMax, biasCorrectionMultiplier));
    return rcpp_result_gen;
END_RCPP
}
// makeEmpDistrList
List makeEmpDistrList(NumericVector X, IntegerMatrix windows, IntegerVector rstSizes, String regridMethod, int maxCore, double fixedMin, double fixedMax, NumericVector biasCorrectionMultiplier);
RcppExport SEXP sourceCpp_58_makeEmpDistrList(SEXP XSEXP, SEXP windowsSEXP, SEXP rstSizesSEXP, SEXP regridMethodSEXP, SEXP maxCoreSEXP, SEXP fixedMinSEXP, SEXP fixedMaxSEXP, SEXP biasCorrectionMultiplierSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type windows(windowsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rstSizes(rstSizesSEXP);
    Rcpp::traits::input_parameter< String >::type regridMethod(regridMethodSEXP);
    Rcpp::traits::input_parameter< int >::type maxCore(maxCoreSEXP);
    Rcpp::traits::input_parameter< double >::type fixedMin(fixedMinSEXP);
    Rcpp::traits::input_parameter< double >::type fixedMax(fixedMaxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type biasCorrectionMultiplier(biasCorrectionMultiplierSEXP);
    rcpp_result_gen = Rcpp::wrap(makeEmpDistrList(X, windows, rstSizes, regridMethod, maxCore, fixedMin, fixedMax, biasCorrectionMultiplier));
    return rcpp_result_gen;
END_RCPP
}
