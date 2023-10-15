// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "h/charlieThreadPool2.hpp"


// [[Rcpp::export]]
SEXP movingAverage(NumericVector x, int windowSize, int speed, 
                   bool returnWindows = true) 
{
  if (speed > x.size() or windowSize > x.size())
    stop("Wrong window size or sliding speed.");
  int rstsize = (x.size() - windowSize + speed) / speed; // Number of windows.
  NumericVector rst(rstsize);
  std::vector<double> C(x.size() + 1);
  std::partial_sum(x.begin(), x.end(), C.begin() + 1);
  C[0] = 0;
  double *cs = C.data() + 1;
  int xsize = x.size();
  const bool returnWin = returnWindows;
  IntegerMatrix wins;
  if (returnWin) wins = IntegerMatrix(2, rstsize);
  for (int i = 0; i < rstsize; ++i)
  {
    int first = i * speed, last = first + windowSize - 1;
    if (last >= xsize) { last = xsize - 1; first = last - windowSize + 1; }
    rst[i] = (cs[last] - cs[first - 1]) / windowSize;
    if (returnWin)
    {
      wins[i * 2] = first + 1;
      wins[i * 2 + 1] = last + 1;
    }
  }
  if (returnWin) return List::create(Named("averages") = rst, 
      Named("windows") = wins);
  return rst;
}


// [[Rcpp::export]]
NumericVector movingAverageSmoothing(NumericVector x, int windowSize, int iterations) 
{
  
  if (windowSize % 2 == 0) stop("Window size even.");
  // const int begin = windowSize / 2, end = rst.size() - begin;
  const int begin = windowSize / 2;
  const int &halfWin = begin;
  std::vector<double> cntr( (x.size() + 2 * halfWin) * 2);
  double *y = cntr.data();
  double *z = y + x.size() + 2 * halfWin;
  const int ysize = x.size() + 2 * halfWin; // &zsize = ysize;
  std::copy(x.begin(), x.end(), y + halfWin);
  for (int i = 0; i < halfWin; ++i)
  {
    z[i] = y[i] = x[0]; 
  }
  for (int i = ysize - 1, iend = ysize - halfWin; i >= iend; --i)
  {
    z[i] = y[i] = *(x.end() - 1);
  }
    
  
  double r = 1.0 / windowSize;
  const int end = ysize - halfWin;
  for (int iter = 0; iter < iterations; ++iter)
  {
    double S = std::accumulate(y, y + windowSize, 0.0);
    for (int i = begin; i < end; ++i)
    {
      z[i] = S * r;
      S += y[i + halfWin] - y[i - halfWin];
    } 
    std::swap(y, z);
  }
  
  
  return NumericVector(y + halfWin, y + end);
} 









/*

 
 
 for (i in 1:30000)
 {
   tmp2 = runif(100)
   tmp = movingAverage(tmp2, sample(100,1), sample(100,1))
   if (max(abs(range(apply(tmp$windows, 2, function(x) 
     mean(tmp2[x[1]:x[2]])) - tmp$averages))) > 1e-10)
   {
     print("wrong!"); break
   }
 }
 
 
 
*/




// [[Rcpp::export]]
DataFrame windowVariances(
    NumericVector x, double windowSizePercentageIncrement = 0.01, 
    int maxCore = 1000)
{
  int sizeIncrement = std::max(
    int(std::round(x.size() * windowSizePercentageIncrement)), 1);
  sizeIncrement = std::min(int(x.size()), sizeIncrement);
  int xsize = x.size();
  struct mm { double s, ss; }; // sum and sum of squares.
  std::vector<mm> C(x.size() + 1);
  C[0].s = C[0].ss = 0;
  mm *cs = C.data() + 1;
  for (int i = 0; i < xsize; ++i)
  {
    cs[i].s = cs[i - 1].s + x[i];
    cs[i].ss = cs[i - 1].ss + x[i] * x[i];
  }
  CharlieThreadPool cp(std::move(maxCore));
  maxCore = cp.maxCore;
  int Njob = xsize / sizeIncrement + 1;
  Njob -= sizeIncrement == 1;
  int grainSize = std::max(1, Njob / (maxCore * maxCore * maxCore));
  NumericVector rst(Njob); 
  IntegerVector winSize(Njob); 
  if (sizeIncrement != 1)
  {
    for (int i = 0, iend = winSize.size(); i < iend; ++i)
      winSize[i] = i * sizeIncrement;
    winSize[0] = 1;
  }
  else
  {
    std::iota(winSize.begin(), winSize.end(), 1);
  }
  double *rst_ = &rst[0];
  int *winSize_ = &winSize[0];
  
  
  auto runner = [winSize_, cs, rst_, xsize](std::size_t i, std::size_t t)->bool
  {
    if (i == 0) return false;
    int windowSize = winSize_[i];
    double S = 0;
    for (int k = windowSize; k <= xsize; ++k)
    {
      int last = k - 1, first = k - windowSize;
      double mean = (cs[last].s - cs[first - 1].s) / windowSize;
      double secondM = (cs[last].ss - cs[first - 1].ss) / windowSize;
      S += secondM - mean * mean;
    }
    rst_[i] = S / (xsize - windowSize + 1);
    return false;
  };
  
  
  cp.parFor(0, Njob, runner, grainSize,
            [](std::size_t t)->bool{ return false; },
            [](std::size_t t)->bool{ return false; });
  
  
  return DataFrame::create(
    Named("windowSize") = winSize, Named("variance") = rst);
}




/*
 
 
 r = exp(runif(1, log(1e-5), log(1)))
 x = runif(1000)
 rst = windowVariances(x, r)
 
 
 truth = sapply(rst$windowSize, function(w)
 {
   a = movingAverage(x, w, 1, F) 
   a2 = movingAverage(x * x, w, 1, F) 
   mean(a2 - a * a)
 })
 
 
 range(truth - rst$variance)
 
 
 
 
 N = 10000
 k = 3000
 x = runif(N)
 ind = sort(sample(N, k))
 x[ind] = sort(x[ind])
 rst = windowVariances(x, 0)
 plot(rst, type = "l")
 
 
 N = 10000
 x = runif(N)
 y = rnorm(N)
 x1 = 2 * x + y; x2 = 3 * x - 2 * y
 x2 = x2[order(x1)]
 rst = windowVariances(x2, 0)
 plot(rst, type = "l")
 
 
 
 
*/















#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// movingAverage
SEXP movingAverage(NumericVector x, int windowSize, int speed, bool returnWindows);
RcppExport SEXP sourceCpp_62_movingAverage(SEXP xSEXP, SEXP windowSizeSEXP, SEXP speedSEXP, SEXP returnWindowsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type windowSize(windowSizeSEXP);
    Rcpp::traits::input_parameter< int >::type speed(speedSEXP);
    Rcpp::traits::input_parameter< bool >::type returnWindows(returnWindowsSEXP);
    rcpp_result_gen = Rcpp::wrap(movingAverage(x, windowSize, speed, returnWindows));
    return rcpp_result_gen;
END_RCPP
}
// movingAverageSmoothing
NumericVector movingAverageSmoothing(NumericVector x, int windowSize, int iterations);
RcppExport SEXP sourceCpp_62_movingAverageSmoothing(SEXP xSEXP, SEXP windowSizeSEXP, SEXP iterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type windowSize(windowSizeSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(movingAverageSmoothing(x, windowSize, iterations));
    return rcpp_result_gen;
END_RCPP
}
// windowVariances
DataFrame windowVariances(NumericVector x, double windowSizePercentageIncrement, int maxCore);
RcppExport SEXP sourceCpp_62_windowVariances(SEXP xSEXP, SEXP windowSizePercentageIncrementSEXP, SEXP maxCoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type windowSizePercentageIncrement(windowSizePercentageIncrementSEXP);
    Rcpp::traits::input_parameter< int >::type maxCore(maxCoreSEXP);
    rcpp_result_gen = Rcpp::wrap(windowVariances(x, windowSizePercentageIncrement, maxCore));
    return rcpp_result_gen;
END_RCPP
}
