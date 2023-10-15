#pragma once

template <typename num>
void enforceNonnegativity(num *x, num *xend)
{
  for (; x < xend; ++x) *x = std::max(0.0, *x);
} 




template <typename ing, typename num1, typename num2,
          typename num3, typename num4,  
          bool addTo, bool zeroFuzzyNegative>
bool LinearRegrid(num1 *x, num2 *xp, ing xsize, 
                  num3 *y, num4 *yp, ing ysize)
{
  yp[0] = 0.0;
  y[0] = y[0] > x[0] ? x[0] : y[0];
  y[ysize - 1] = y[ysize - 1] < x[xsize - 1] ? x[xsize - 1] : y[ysize - 1];
  double tmp = 0, interval = 0;
  ing j = 0;
  for (ing i = 1; i < ysize; ++i)
  {
    if (!addTo) yp[i] = 0; 
    interval = y[i] - (double)y[i - 1];
    for (; j < xsize and x[j] <= y[i]; ++j)
    {
      tmp = xp[j] / interval * ((double)y[i] - x[j]);
      yp[i - 1] += tmp;
      yp[i] += xp[j] - tmp;
      if (zeroFuzzyNegative) yp[i] = std::max(yp[i], 0.0);
    }
  }
  return false; 
}






template <
  typename ing, typename num1, typename num2,
  typename num3, typename num4, bool addTo, 
  bool zeroFuzzyNegative>
bool FourPointRegrid(
    num1 *inLoss,  num2 *inProb,  ing inSize, 
    num3 *outLoss, num4 *outProb, ing outSize,
    std::vector<double> &safe)
{
  
  
  bool unableToWork = false;
  if (addTo) safe.assign(outProb, outProb + outSize);
  
  
  outLoss[0] = outLoss[0] > inLoss[0] ? inLoss[0] : outLoss[0];
  outLoss[outSize - 1] = outLoss[outSize - 1] < inLoss[inSize - 1] ? 
  inLoss[inSize - 1] : outLoss[outSize - 1];
  
  
  const double EPS100K = 1e-12;
  
  
  if (outSize < 5 or inSize < 5) unableToWork = true;
  else
  {
    if (!addTo) std::fill(outProb, outProb + outSize, 0.0);
    ing n = 0;
    
    
    double x0 = outLoss[0];
    double x3 = outLoss[outSize - 1];
    double p0 = 0;
    double p1 = 0;
    double p2 = 0;
    double p3 = 0;
    double p = 0;
    
    
    double x = 0;
    double x1 = 0;
    double x2 = 0;
    
    
    for (ing i = 0; i < outSize - 1; ) 
    {
      if (std::abs(inLoss[n] - outLoss[i]) < EPS100K) 
      {
        outProb[i] += (double)inProb[n];
        ++n;
        if (n >= inSize) break;
      } 
      else 
      {
        x1 = outLoss[i];
        x2 = outLoss[i + 1];
        
        
        while (n < inSize and (x2 - inLoss[n]) > EPS100K)
        {
          x = inLoss[n];
          p = inProb[n];
          double tmp = p * (x3 - x) * (x0 - x) /
            ((x3 - x1) * (x0 - x1) * (x2 - x) +
              (x3 - x2) * (x0 - x2) * (x - x1));
          p1 = tmp * (x2 - x);
          p2 = tmp * (x - x1);
          p0 = (p * (x3 - x) - (x3 - x1) * p1 - (x3 - x2) * p2) / (x3 - x0);
          p3 = p - p0 - p1 - p2;
          ++n;
          outProb[0] += (double)p0;
          outProb[i]     += std::max(0.0, (double)p1); // outProb[i] += (double)p1;
          outProb[i + 1] += std::max(0.0, (double)p2); // outProb[i + 1] += (double)p2;
          outProb[outSize - 1] += (double)p3;
        }
        ++i;
      }
    }
    
    
    if (std::abs(outLoss[outSize - 1] - inLoss[inSize - 1]) < EPS100K)
      outProb[outSize - 1] += inProb[inSize - 1];
    
    
    if (outProb[0] >= 0 and outProb[outSize - 1] >= 0) 
    {
      if (zeroFuzzyNegative) enforceNonnegativity(outProb, outProb + outSize);
      return 0;
    }
    
    
    ing i = 0;
    ing j = outSize - 1;
    
    
    while (true) 
    {
      if (outProb[i] < 0) 
      {
        x = outLoss[i];
        x1 = outLoss[i + 1];
        x2 = outLoss[i + 2];
        x3 = outLoss[j];
        p = outProb[i];
        p1 = (x - x3) * (x - x2) / ((x1 - x3) * (x1 - x2)) * p;
        p2 = (x - x1) * (x - x3) / ((x2 - x3) * (x2 - x1)) * p;
        p3 = p - p1 - p2;
        outProb[i + 1] += (double)p1;
        outProb[i + 2] += (double)p2;
        outProb[j] += (double)p3;
        outProb[i] = 0;
        ++i;
      }
      
      
      if (outProb[i] >= 0 and outProb[j] >= 0) break;
      
      
      if (j - i < 3)
      {
        unableToWork = true;
        break;
      }
      
      
      if (outProb[j] < 0) 
      {
        x = outLoss[j];
        x1 = outLoss[j - 1];
        x2 = outLoss[j - 2];
        x3 = outLoss[i];
        p = outProb[j];  
        p1 = (x - x3) * (x - x2) / ((x1 - x3) * (x1 - x2)) * p;
        p2 = (x - x1) * (x - x3) / ((x2 - x3) * (x2 - x1)) * p;
        p3 = p - p1 - p2;
        outProb[j - 1] += (double)p1;
        outProb[j - 2] += (double)p2;
        outProb[i] += (double)p3;
        outProb[j] = 0;
        --j;
      }
      
      
      if (outProb[i] >= 0 and outProb[j] >= 0) break;
      
      
      if (j - i < 3)  
      {
        unableToWork = true;
        break;
      }
    }
  }
  
  
  if (unableToWork)
  {
    if (addTo) std::copy(safe.begin(), safe.end(), outProb);
    LinearRegrid <ing, num1, num2, num3, num4, addTo, zeroFuzzyNegative> (
        inLoss, inProb, inSize, outLoss, outProb, outSize);
  }
  else
  {
    if (zeroFuzzyNegative) enforceNonnegativity(outProb, outProb + outSize);
  }
  
  
  return unableToWork;
}






template <
  typename ing, typename num1, typename num2, 
  typename num3, typename num4, 
  bool addTo, bool zeroFuzzyNegative>
bool LMM(
    num1 *x,    num2 *p,    ing siz, 
    num3 *xnew, num4 *pnew, ing sizeNew,
    std::vector<double> &safe)
{
  
  // if (xnew[sizeNew - 1] < x[siz - 1]) xnew[sizeNew - 1] = x[siz - 1];
  // if (xnew[0] > x[0]) xnew[0] = x[0];
  xnew[0] = xnew[0] > x[0] ? x[0] : xnew[0];
  xnew[sizeNew - 1] = 
    xnew[sizeNew - 1] < x[siz - 1] ? x[siz - 1] : xnew[sizeNew - 1];
  
  
  num3 *left = &xnew[0], *xnewEnd = xnew + sizeNew;
  num4 *leftP = &pnew[0];
  
  
  ing i = 0, iend = siz;
  double x1_x0 = left[1] - left[0], x2_x1 = left[2] - left[1];
  double x2_x0 = x2_x1 + x1_x0, leftHeadTailMid = (left[0] + left[3]) / 2;
  
  
  if (addTo) safe.assign(pnew, pnew + sizeNew); 
  else std::fill(pnew, pnew + sizeNew, 0.0);
  
  
  bool notWorking = false;
  
  
  if (siz < 4 or sizeNew < 4) notWorking = true;
  else 
  {
    for (; i < iend; ++i) 
    {
      if (x[i] > leftHeadTailMid)  
      {
        
        
        if (leftP[0] < 0) 
        {
          notWorking = 1;
          break;
        }
        ++left;
        
        
        if (left + 2 >= xnewEnd) 
        {
          --left;
          break;
        }
        
        
        ++leftP;
        
        
        leftHeadTailMid = (left[0] + left[3]) / 2;
        x1_x0 = x2_x1;
        x2_x1 = left[2] - left[1];
        x2_x0 = x2_x1 + x1_x0;
      }
      
      
      leftP[0] += (x[i] - left[1]) * (x[i] - left[2]) / (x1_x0 * x2_x0) * p[i];
      leftP[1] += (x[i] - left[0]) * (x[i] - left[2]) / (-x1_x0 * x2_x1) * p[i];
      leftP[2] += (x[i] - left[0]) * (x[i] - left[1]) / (x2_x0 * x2_x1) * p[i];
    }
  }
  
  
  if (!notWorking and i != iend) 
  {
    for (; i < iend; ++i) 
    {
      leftP[0] += (x[i] - left[1]) * (x[i] - left[2]) / (x1_x0 * x2_x0) * p[i];
      leftP[1] += (x[i] - left[0]) * (x[i] - left[2]) / (-x1_x0 * x2_x1) * p[i];
      leftP[2] += (x[i] - left[0]) * (x[i] - left[1]) / (x2_x0 * x2_x1) * p[i];
    }
  }
  
  
  if (!notWorking && leftP[0] >= 0 && leftP[1] >= 0 && leftP[2] >= 0) 
  {
    if (zeroFuzzyNegative) enforceNonnegativity(pnew, pnew + sizeNew);
    return false;
  }
  
  
  // probability container polluted, copy the safe and start linear regridding
  if (addTo) std::copy(safe.begin(), safe.end(), pnew);
  
  
  LinearRegrid<ing, num1, num2, num3, num4, addTo, zeroFuzzyNegative> (
      x, p, siz, xnew, pnew, sizeNew);
  
  
  return true;
} 






// method == 0: linear regrid.
// method == 1: local moment matching.
// method == 2: 4-point regrid.
struct Regrid
{
  std::vector<double> safe;
  
  
  template <int method, typename ing, 
            typename num1, typename num2,
            typename num3, typename num4, 
            bool addTo, bool zeroFuzzyNegative>
  bool operator()(
      num1 *x,    num2 *p,    ing siz, 
      num3 *xnew, num4 *pnew, ing sizeNew
  )
  {
    
    if (method == 0) return LinearRegrid <ing, 
      num1, num2, num3, num4, addTo, zeroFuzzyNegative>(
          x, p, siz, xnew, pnew, sizeNew
      );
    
    
    if (method == 1) return LMM <ing,
      num1, num2, num3, num4, addTo, zeroFuzzyNegative>(
          x, p, siz, xnew, pnew, sizeNew, safe
      );
    
    
    if (method == 2) return FourPointRegrid <ing,
      num1, num2, num3, num4, addTo, zeroFuzzyNegative>(
          x, p, siz, xnew, pnew, sizeNew, safe
      );
    
    
    // static_assert(method >= 0 and method <= 2);
    
    
    return false;
  }
};




































