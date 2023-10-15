// #pragma once


#include <Rcpp.h>
using namespace Rcpp;

unsigned Bsrch(double value, NumericVector::iterator low, NumericVector::iterator up)
{
  NumericVector::iterator st = low, mid;//store the start of the sequence
  --up;
  if(value <= *low)return 0;
  if(value >* up)return up - st;
  while(true)
  {
    mid = low + unsigned((up - low) / 2);
    if(value <= *mid) up = mid;
    else if(value > *(mid + 1)) low = mid;
    else return up - st;
  }
}


// [[Rcpp::export]]
std::vector<double> mapTo(NumericVector rv, DataFrame cdf)
{
  NumericVector val=cdf[0], P=cdf[1];
  NumericVector::iterator ri=rv.begin(), valst=val.begin();
  std::vector<double>rst(rv.size());
  std::vector<double>::iterator rsti=rst.begin();
  for(;ri!=rv.end();++ri,++rsti)
    *rsti=*(valst+Bsrch(*ri,P.begin(),P.end()));
  return rst;
}


// [[Rcpp::export]]
NumericVector LHS(DataFrame pdf, int N)
{
  NumericVector val = pdf[0], P = pdf[1];
  int dsiz = val.size();
  std::vector<double> cumP(dsiz);
  std::partial_sum(P.begin(), P.end(), cumP.begin());
  cumP.back() = 1;
  NumericVector r = Language("runif", N).eval();
  double gap = 1.0 / N;
  for(int i = 0, iend = r.size(); i < iend; ++i)
  {
    r[i] = (r[i] + i) * gap;
  }
  std::vector<double> rst(N);
  int j = 0;
  for(int i = 0; i < dsiz; ++i)
  {
    while(r[j] <= cumP[i] and j < N)
    {
      rst[j] = val[i];
      ++j;
    }
  }
  IntegerVector shuf = Language("sample", N).eval();
  NumericVector rstR(N);
  for(int i = 0; i < N; ++i)
  {
    rstR[i] = rst[shuf[i] - 1];
  }
  return rstR;
}

















