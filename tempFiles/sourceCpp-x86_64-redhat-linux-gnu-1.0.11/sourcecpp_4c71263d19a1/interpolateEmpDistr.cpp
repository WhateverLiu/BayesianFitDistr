// [Rcpp::depends(cpp17)]
#include <Rcpp.h>
using namespace Rcpp;
#include "h/mixRegrid.h"
#include "h/charlieThreadPool2.hpp"


template <int regridMethod>
List findEmpDistrGivenMDR(List X, NumericVector MDR, 
                          NumericVector MDRwanted, 
                          IntegerVector sizeWanted,
                          int maxCore)
{
  if (X.size() != MDR.size()) stop("X.size() != MDR.size().");
  
  
  if (sizeWanted.size() != 1 and sizeWanted.size() != MDRwanted.size()) 
    stop("sizeWanted.size() != MDRwanted.size().");
  
  
  int NdistrWanted = MDRwanted.size();
  int givenDistrListSize = X.size();
  
  
  for (int i = 1, iend = MDR.size(); i < iend; ++i)
  {
    if (MDR[i] < MDR[i - 1]) stop(
        "MDR is not sorted. The PMF list X should also have been sorted along with MDR.");
  } 
  
  
  for (int i = 0, iend = MDRwanted.size(); i < iend; ++i)
  {
    if (MDRwanted[i] < MDR[0] or MDRwanted[i] > *(MDR.end() - 1))
      stop("MDRwanted are not fully covered by MDR.");
  }
  
  
  List rst(NdistrWanted);
  struct distr { int size; double *val, *P; };
  std::vector<distr> rstDistrCntr(NdistrWanted + X.size());
  distr *rstDistrList = rstDistrCntr.data();
  distr *givenDistrList = rstDistrList + NdistrWanted;
  
  
  int sz = sizeWanted.size();
  for (int i = 0; i < NdistrWanted; ++i)
  {
    NumericVector val(sizeWanted[i % sz]), P(sizeWanted[i % sz]);
    rstDistrList[i].val = &val[0];
    rstDistrList[i].P = &P[0];
    rstDistrList[i].size = sizeWanted[i % sz];
    rst[i] = List::create(Named("val") = val, Named("P") = P);
  }
  
  
  for (int i = 0; i < givenDistrListSize; ++i)
  {
    List tmp = X[i];
    NumericVector val = tmp[0], P = tmp[1];
    givenDistrList[i].val = &val[0];
    givenDistrList[i].P = &P[0];
    givenDistrList[i].size = val.size();
  }
  
  
  double *mdr = &MDR[0], *mdrWanted = &MDRwanted[0];
  
  
  CharlieThreadPool cp(std::move(maxCore));
  // maxCore = cp.maxCore;
  std::vector<MixRegrid> MR(maxCore);
  
  
  // Rcout << "1.1\n";
  
  
  auto f = [&](std::size_t i, std::size_t t)->bool
  {
    MixRegrid &mr = MR[t];
    auto it = std::lower_bound(mdr, mdr + givenDistrListSize, mdrWanted[i]);
    double left = *(it - 1), right = *it;
    double w = (mdrWanted[i] - left) / (right - left);
    int rightInd = it - mdr;
    int leftInd = rightInd - 1;
    
    
    mr.operator()<regridMethod, int, double, double, double, double, double, double> (
        givenDistrList[leftInd].val ,  givenDistrList[leftInd].P, givenDistrList[leftInd].size, 
        givenDistrList[rightInd].val, givenDistrList[rightInd].P, givenDistrList[rightInd].size, 
        rstDistrList[i].val, rstDistrList[i].P, rstDistrList[i].size, w);
    
    
    return false;
  };
  
  
  int grainSize = NdistrWanted / (maxCore * maxCore) + 1;
  cp.parFor(0, NdistrWanted, f, grainSize, 
            [](std::size_t t)->bool { return false; },
            [](std::size_t t)->bool { return false; });
  
  
  return rst;
}




// [[Rcpp::export]]
List findEmpDistrGivenMDR(List X, NumericVector MDR, 
                          NumericVector MDRwanted, 
                          IntegerVector sizeWanted,
                          int maxCore = 1000,
                          String regridMethod = "r4")
{
  if (regridMethod == "lr")
    return findEmpDistrGivenMDR<0> (X, MDR, MDRwanted, sizeWanted, maxCore);
  
  
  if (regridMethod == "lmm")
    return findEmpDistrGivenMDR<1> (X, MDR, MDRwanted, sizeWanted, maxCore);
  
  
  if (regridMethod == "r4")
    return findEmpDistrGivenMDR<2> (X, MDR, MDRwanted, sizeWanted, maxCore);
  
  
  stop("Regrid method not implemented.");
  
  
  return List::create();
}









































#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// findEmpDistrGivenMDR
List findEmpDistrGivenMDR(List X, NumericVector MDR, NumericVector MDRwanted, IntegerVector sizeWanted, int maxCore, String regridMethod);
RcppExport SEXP sourceCpp_54_findEmpDistrGivenMDR(SEXP XSEXP, SEXP MDRSEXP, SEXP MDRwantedSEXP, SEXP sizeWantedSEXP, SEXP maxCoreSEXP, SEXP regridMethodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MDR(MDRSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MDRwanted(MDRwantedSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sizeWanted(sizeWantedSEXP);
    Rcpp::traits::input_parameter< int >::type maxCore(maxCoreSEXP);
    Rcpp::traits::input_parameter< String >::type regridMethod(regridMethodSEXP);
    rcpp_result_gen = Rcpp::wrap(findEmpDistrGivenMDR(X, MDR, MDRwanted, sizeWanted, maxCore, regridMethod));
    return rcpp_result_gen;
END_RCPP
}
