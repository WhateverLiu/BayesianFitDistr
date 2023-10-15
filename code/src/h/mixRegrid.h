




// [x, x + xsize) stores unsorted samples of distribution X.
// [y, y + ysize) and [yp, yp + ysize) store PMF of Y.
// Goal: make a mixture PMF of X and Y. The weight on Y is yweight.
struct MixRegrid
{
  Regrid Regrider;
  std::vector<double> cntr;
  
  
  template <int method, typename ing, typename num1, typename num2,
            typename num3, typename num4, typename num5, typename num6>  
  bool operator()(
      num1 *x,   num2 *xp,   ing xsize, // Sampled sequence for x.
      num3 *y,   num4 *yp,   ing ysize,
      num5 *rst, num6 *rstP, ing rstSize, // Result PMF is written to rst and rstP.
      double yweight)
  {
    
    bool b = false;
    if (yweight >= 1)
    {
      makeRegularGrid(y[0], y[ysize - 1], rst, rstSize); 
      b = Regrider.operator()<
        method, ing, num3, num4, num5, num6, false, true> (
        y, yp, ysize, rst, rstP, rstSize);
    }
    else if (yweight <= 0)
    {
      makeRegularGrid(x[0], x[xsize - 1], rst, rstSize); 
      b = Regrider.operator()<
        method, ing, num1, num2, num5, num6, false, true> (
            x, xp, xsize, rst, rstP, rstSize);
    }
    else 
    {
      ing mergeSize = xsize + ysize;
      cntr.resize( mergeSize * 2 );
      double *merged = cntr.data();
      double *mergedP = merged + mergeSize;
      pmfMerge(x, xp, xsize, y, yp, ysize, merged, mergedP, yweight); 
      pmfCollapse(merged, mergedP, mergeSize);
      makeRegularGrid(merged[0], merged[mergeSize - 1], rst, rstSize);
      b = Regrider.operator()<method, ing, double, double, num5, num6, false, true> (
        merged, mergedP, mergeSize, rst, rstP, rstSize);
    }
    
    
    return b;
  }
  
  
  
  
  
  
  
  
};
















