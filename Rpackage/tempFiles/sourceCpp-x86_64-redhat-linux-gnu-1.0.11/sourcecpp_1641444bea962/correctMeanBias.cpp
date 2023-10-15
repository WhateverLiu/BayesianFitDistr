// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "h/charlieThreadPool2.hpp"
#include "h/upscaleBoundedMean.h"


/*
struct MeanCorrection
{
  std::vector<double> cntr;
  // [x, x + size) has been sorted. upperBound > targetMean.
  double operator() (const double *x_, const int size, const double targetMean,
                  const double upperBound, const double eps)
  {
    // Compute the ratio via sorting. Probably will be slower than the plain iterative algo.
    if (false)
    {
      cntr.resize(size + 1 + size);
      double *x = cntr.begin() + size + 1;
      std::copy(x_, size, x);
      std::sort(x, x + size);
      cntr[0] = 0;
      double *cs = cntr.data() + 1;
      std::partial_sum(x, x + size, cs);
      double r = targetMean / cs[size - 1];
      while (true)
      {
        // How many elements would be >= upperBound.
        int howManyExceeded = x + size - std::lower_bound(
          x, x + size, upperBound, [r](const double &u, const double &val)->bool
          {
            return u * r < val;
          });
        int howManyBelow = size - howManyExceeded;
        double meanWouldBe = (cs[howManyBelow] + howManyExceeded) / size;
        if (std::abs(meanWouldBe - targetMean) < eps) break;
        r *= targetMean / meanWouldBe;
      }
      return r;
    }
    
    
    double *x = x_, *xend = x + size;
    double r = 1;
    while (true)
    {
      double mean = 0;
      for (auto it = x; it < xend; ++it) mean += std::min(*it * r, upperBound);
      mean /= size;
      if (std::abs(mean - targetMean) < eps) break;
      r *= targetMean / mean;
    }
    return r;
  }
};


double meanCorrection(const double *x_, const int size, 
                      const double targetMean, 
                      const double upperBound, const double eps)
{
  double *x = x_, *xend = x + size;
  double r = 1;
  while (true)
  {
    double mean = 0;
    for (auto it = x; it < xend; ++it) mean += std::min(*it * r, upperBound);
    mean /= size;
    if (std::abs(mean - targetMean) < eps) break;
    r *= targetMean / mean;
  }
  return r;
}


// windows are 1-based inclusive intervals.
// [[Rcpp::export]]
NumericVector meanCorrectionScaler(
    NumericVector mdrs, NumericVector cdrs, IntegerMatrix windows,
    double lim = 1, double eps = 1e-10, int maxCore = 1000) 
{
  if (mdrs.size() != cdrs.size()) stop("mdrs.size() != cdrs.size().");
  int *wins = &windows[0];
  int Nwins = windows.ncol();
  charlieThreadPool cp(maxCore);
  maxCore = cp.maxCore;
  NumericVector rst(Nwins);
  auto f = [&](std::size_t i, std::size_t t)->bool
  {
    int begin = wins[i * 2] - 1, end = wins[i * 2 + 1], size = end - begin;
    double targetMean = std::accumulate(
      mdrs.begin() + begin, mdrs.begin() + end, 0.0) / size;
    if (targetMean > lim) rst[i] = 1e300;
    else rst[i] = meanCorrection(cdrs.begin(), size, targetMean, lim, eps);
    return false;
  };
  int grainSize = Nwins / (maxCore * maxCore * maxCore) + 1;
  cp.parFor(0, Nwins, f, grainSize, [](std::size_t t)->bool { return false; },
            [](std::size_t t)->bool { return false; });
  for (auto it = rst.begin(); it < rst.end(); ++it)
    if (*it >= 1e300) stop("Target mean greater than upper bound.");
  return rst;
}
*/




// [[Rcpp::export]]
List correctMeanBias(List distlist, NumericVector mdrs, double lim = 1.0, 
                     double eps = 1e-10, int maxCore = 1000)
{
  struct dist {int size; double *val, *P;};
  int Ndist = distlist.size();
  if (mdrs.size() != Ndist) stop("distlist.size() != mdrs.size()");
  std::vector<dist> D(Ndist);
  for (int i = 0; i < Ndist; ++i)
  {
    List d = distlist[i];
    NumericVector val = d[0], p = d[1];
    D[i].size = val.size();
    D[i].val = &val[0];
    D[i].P = &p[0];
  }
  
  
  /*
  auto findScaler = [](
    double *val, double *P, int size, double targetMean, 
    double upperBound, double eps)->double
  {
    double r = 1;
    while (true)
    {
      double mean = 0;
      for (int i = 0; i < size; ++i)
        mean += std::min(val[i] * r, upperBound) * P[i];
      if (std::abs(mean - targetMean) < eps) break;
      r *= targetMean / mean;
    }
    return r;
  };
  */
  
  
  CharlieThreadPool cp(std::move(maxCore));
  // maxCore = cp.maxCore;
  // std::vector<UpscaleBoundedMean> UBM(maxCore);
  
  
  int grainSize = Ndist / (maxCore * maxCore * maxCore) + 1;
  std::vector<double> multipliers(Ndist);
  cp.parFor(0, Ndist, [&](std::size_t i, std::size_t t)->bool
  {
    if (mdrs[i] > lim) multipliers[i] = 1e300;
    else
    {
      int iter = 0;
      // multipliers[i] = UBM[t](D[i].val, D[i].P, D[i].size, mdrs[i], lim, iter);
      multipliers[i] = upscaleBoundedMean(
        D[i].val, D[i].P, D[i].size, mdrs[i], lim, eps, iter, 100);
    }
    // multipliers[i] = findScaler(
    //   D[i].val, D[i].P, D[i].size, mdrs[i], lim, eps);
    return false;
  }, grainSize, [](std::size_t t)->bool {return false;}, 
  [](std::size_t t)->bool {return false;});
  
  
  for (auto it = multipliers.begin(); it < multipliers.end(); ++it)
  {
    if (*it >= 1e300) 
    {
      warning("Target mean >= upperBound. Skip correction and set the multiplier for scaling PMF support to 1.");
      *it = 1.0;
    }
  }
    
  
  List rst(Ndist);
  for (int i = 0; i < Ndist; ++i)
  {
    List d = distlist[i];
    NumericVector val = d[0], P = d[1];
    int j = val.size() - 1;
    double mr = multipliers[i];
    for (; j >= 0; --j)
    {
      if (val[j] * mr < lim) break;
    }
    if (j >= val.size() - 1) // PMF needs to be maintained or scaled down.
    {
      if (std::abs(mr - 1.0) < eps) rst[i] = d; // If scaler is effectively 1.0. PMF should be maintained.
      else
      {
        NumericVector valnew(val.size());
        for (int k = 0, kend = val.size(); k < kend; ++k) 
          valnew[k] = val[k] * mr;
        rst[i] = List::create(Named("val") = valnew, Named("P") = P);   
      }
    }
    else // PMF needs to be scaled up.
    {
      // j is the index of the last value that would be < lim after adjustment.
      int sizenew = j + 2;
      NumericVector valnew(sizenew);
      for (int k = 0, kend = sizenew - 1; k < kend; ++k) valnew[k] = val[k] * mr;
      valnew[sizenew - 1] = lim;
      NumericVector Pnew(sizenew);
      double psum = 0;
      for (int k = 0, kend = sizenew - 1; k < kend; ++k)
      {
        psum += P[k];
        Pnew[k] = P[k];
      }
      Pnew[sizenew - 1] = 1.0 - psum;
      rst[i] = List::create(Named("val") = valnew, Named("P") = Pnew);   
    }
    
    
  }
  
  
  return rst;
}


























#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// correctMeanBias
List correctMeanBias(List distlist, NumericVector mdrs, double lim, double eps, int maxCore);
RcppExport SEXP sourceCpp_5_correctMeanBias(SEXP distlistSEXP, SEXP mdrsSEXP, SEXP limSEXP, SEXP epsSEXP, SEXP maxCoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type distlist(distlistSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mdrs(mdrsSEXP);
    Rcpp::traits::input_parameter< double >::type lim(limSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxCore(maxCoreSEXP);
    rcpp_result_gen = Rcpp::wrap(correctMeanBias(distlist, mdrs, lim, eps, maxCore));
    return rcpp_result_gen;
END_RCPP
}
