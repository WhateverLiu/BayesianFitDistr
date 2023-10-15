#pragma once
#include "processPmfs.h"
#include "regrid.h"
// [x, x + xsize) stores unsorted samples of distribution X.
// [y, y + ysize) and [yp, yp + ysize) store PMF of Y.
// Goal: make a mixture PMF of X and Y. The weight on Y is yweight.
struct MergeRegrid
{
  Regrid Regrider;
  std::vector<double> cntr;
  
  
  template <int method, typename ing, typename num1, typename num2,
            typename num3, typename num4, typename num5, typename num6>  
  bool operator()(
      num1 *x,   ing xsize, // Sampled sequence for x.
      num2 *y,   num3 *yp,   ing ysize,
      num4 *rst, num5 *rstP, ing rstSize, // Result PMF is written to rst and rstP.
      double yweight, double fixedMin, double fixedMax,
      num6 biasCorrectionMultiplier)
  {
    
    bool b = false;
    num6 bm = biasCorrectionMultiplier;
    if (yweight >= 1)
    {
      
      double minval = fixedMin == 1e300 ? y[0] : fixedMin;
      double maxval = fixedMax == 1e300 ? y[ysize - 1] : fixedMax;
      makeRegularGrid(minval, maxval, rst, rstSize);
      b = Regrider.operator()<
        method, ing, num2, num3, num4, num5, false, true> (
        y, yp, ysize, rst, rstP, rstSize);
      
    }
    else if (yweight <= 0)
    {
      cntr.resize(xsize * 2);
      double *u = cntr.data(), *up = u + xsize;
      // std::copy(x, x + xsize, u);
      for (ing i = 0; i < xsize; ++i) u[i] = x[i] * bm;
      ing usize = xsize;
      std::sort(u, u + usize);
      std::fill(up, up + usize, 1.0 / usize);
      pmfCollapse(u, up, usize);
      
      
      double minval = fixedMin == 1e300 ? u[0] : fixedMin;
      double maxval = fixedMax == 1e300 ? u[usize - 1] : fixedMax;
      makeRegularGrid(minval, maxval, rst, rstSize);
      
      
      b = Regrider.operator()<
        method, ing, double, double, num4, num5, false, true> (
            u, up, usize, rst, rstP, rstSize);
    }
    else 
    {
      
      
      ing mergeSize = xsize + ysize;
      cntr.resize( xsize + xsize + mergeSize * 2 );
      double *u = cntr.data(), *up = u + xsize;  
      ing usize = xsize;
      // std::copy(x, x + xsize, u);
      for (ing i = 0; i < xsize; ++i) u[i] = x[i] * bm;
      std::sort(u, u + usize);
      std::fill(up, up + usize, 1.0 / usize);
      double *merged = up + usize, *mergedP = merged + mergeSize;
      pmfMerge(u, up, usize, y, yp, ysize, merged, mergedP, yweight);
      pmfCollapse(merged, mergedP, mergeSize);
      
      
      double minval = fixedMin == 1e300 ? merged[0] : fixedMin;
      double maxval = fixedMax == 1e300 ? merged[mergeSize - 1] : fixedMax;
      makeRegularGrid(minval, maxval, rst, rstSize);
      
      
      b = Regrider.operator()<method, ing, double, double, num4, num5, false, true> (
        merged, mergedP, mergeSize, rst, rstP, rstSize);
    }
    
    
    return b;
  }
  
  
};
















