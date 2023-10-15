#pragma once
#include "genericLongestMonoSubseq.hpp"
#define vec std::vector


// =============================================================================
// Given y[size] and SORTED x[size], find the longest monotonic 
//   subsequence of y[] and the associated subsequence of x. Linearly 
//   interpolate a function y = f(x).
// The class should overload operator() for query.
//
// monoType == 0: increasing subsequence.
// monoType == 1: nondecreasing subsequence.
// monoType == 2: decreasing subsequence.
// monoType == 3: nonincreasing subsequence.
// =============================================================================
template <typename ing, typename numX, typename numY>
struct MonoLinearInterpo
{
  vec<numX> xsubseq;
  vec<numY> ysubseq;
  vec<ing> indx;
  LongestMonoSubseq<ing, numY> lms;
 
 
  template<int monoType>
  void reset(numX *x, numX *xend, numY *y)
  {
    indx.resize(xend - x);
    auto indxEnd = lms.template operator()<monoType>(
      y, y + (xend - x), indx.data());
    indx.resize(indxEnd - indx.data());
    xsubseq.resize(indx.size());
    ysubseq.resize(indx.size());
    for (int64_t i = 0, iend = indx.size(); i < iend; ++i)
    {
      xsubseq[i] = x[indx[i]];
      ysubseq[i] = y[indx[i]];
    }
  } 
  
  
  numY operator()(numX &&x)
  {
    auto it = std::lower_bound(xsubseq.begin(), xsubseq.end(), x);
    if (it >= xsubseq.end()) return ysubseq.back();
    if (it <= xsubseq.begin()) return ysubseq.front();
    auto xlow = *(it - 1), xhigh = *it;
    double w = (x - xlow) / (xhigh - xlow);
    auto yi = it - xsubseq.begin();
    return (1 - w) * ysubseq[yi - 1] + w * ysubseq[yi];
  }
  
  
  numY operator()(numX &x) { return (*this)(std::move(x));  }
  
  
  void operator()(numX *x, numX *xend, numY *y)
  {
    for (auto i = x-x, iend = xend - x; i < iend; ++i)
      y[i] =
        (*this)(x[i]);
  }
};



#undef vec







