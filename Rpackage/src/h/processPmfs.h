#pragma once


// Remove the duplicates. x is sorted.
template <typename ing, typename num1, typename num2>
void pmfCollapse(num1 *x, num2 *p, ing &size)
{
  ing i = 0, j = 1;
  while (j < size)
  {
    if (x[i] == x[j]) p[i] += p[j];
    else
    {
      i += 1; 
      x[i] = x[j];
      p[i] = p[j];
    }
    j += 1;
  }
  size = i + 1;
}


template <typename ing, 
          typename num1, typename num2, 
          typename num3, typename num4,
          typename num5, typename num6>
void pmfMerge(num1 *x, num2 *xp, ing xsize,
              num3 *y, num4 *yp, ing ysize,
              num5 *merged, num6 *mergedP,
              double yweight)
{
  ing i = 0, j = 0, k = 0;
  while (true)
  {
    if (j >= ysize)
    {
      for (; i < xsize; ++i, ++k) 
      { merged[k] = x[i]; mergedP[k] = xp[i] * (1 - yweight); }
      break;
    }
    if (i >= xsize)
    {
      for (; j < ysize; ++j, ++k) 
      { merged[k] = y[j]; mergedP[k] = yp[j] * yweight; }
      break;
    }
    
    
    if (x[i] < y[j]) // or if j >= ysize and i < usize
    {
      merged [k] = x[i]; 
      mergedP[k] = xp[i] * (1.0 - yweight);
      i += 1;
    }
    else
    {
      merged [k] = y[j]; 
      mergedP[k] = yp[j] * yweight;
      j += 1;
    }
    k += 1;
  }
}


template <typename ing, typename num>
void makeRegularGrid(double minval, double maxval, num *rst, ing size)
{
  double delta = (maxval - minval) / (size - 1);
  rst[0] = minval;
  for (ing i = 1, iend = size - 1; i < iend; ++i)
    rst[i] = minval + i * delta;
  rst[size - 1] = maxval;
}



































