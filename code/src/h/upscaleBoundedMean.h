

/*
 inline double scaledBoundedMean(double r, double *val, double *P, 
 int size, double lim)
 {
 double S= 0;
 for (int i = 0; i < size; ++i)
 S += std::min(val[i] * r, lim) * P[i];
 return S;
 }
 inline double upscaleBoundedMean(
 double *val, double *P, int size, double targetMean, 
 double valUB, double eps, int &iter, int maxIter)
 {
 double mean = 0;
 for (int i = 0; i < size; ++i) mean += val[i] * P[i];
 if (mean > targetMean) return targetMean / mean;
 
 
 r = targetMean / mean;
 iter = 0;
 while (iter < maxIter)
 {
 iter += 1;
 mean = scaledBoundedMean(r, val, P, size, valUB);
 if (std::abs(mean - targetMean) < eps) break;
 r *= targetMean / mean;
 }
 return r;
 };  
 */




inline double upscaleBoundedMean(
    double *val, double *P, int size, double targetMean, 
    double valUB, double eps, int &iter, int maxIter)
{
  // double mean = 0;
  // for (int i = 0; i < size; ++i) mean += val[i] * P[i];
  double mean = std::inner_product(val, val + size, P, 0.0);
  if (mean > targetMean) return targetMean / mean;
  
  
  auto f = [&](double r)->std::pair<double,double>
  {
    double m = 0, deriv = 0;
    for (int i = 0; i < size; ++i) 
    {
      double v = val[i] * r;
      if (v < valUB) 
      {
        // double v = val[i] * P[i];
        m += v * P[i];
        deriv += val[i] * P[i];
      }
      else m += valUB * P[i];
    }
    return std::pair<double, double>(m - targetMean, deriv);
  };
  
  
  return newtonRoot(iter, 1.0001, f, eps, maxIter);
}; 







// Redundant.
struct UpscaleBoundedMean
{
  std::vector<double> cntr;
  

  template<typename ing, typename num>  
  double operator()(num *val, num *P, ing size, 
                  num targetMean, num valUB, ing &iter)
  {
    double mean = 0;
    for (ing i = 0; i < size; ++i) mean += val[i] * (double)P[i];
    if (mean > targetMean) return targetMean / mean;
    
    
    cntr.resize((size + 1) * 2);
    double *psum = cntr.data();
    psum[0] = 0;
    psum += 1;
    double *msum = psum + size;
    msum[0] = 0;
    msum += 1;
    for (ing i = 0; i < size; ++i)
    {
      psum[i] = psum[i - 1] + P[i];
      msum[i] = msum[i - 1] + P[i] * val[i];
    }
    
    
    iter = 0;
    auto f = [&](const num &vi, double targetM)->bool
    {
      iter += 1;
      ing i = &vi - val; // When vi is the last value point.
      double lastP = 1.0 - psum[i - 1];
      double scaler = valUB / vi;
      double partialMean = msum[i - 1] * scaler;
      return partialMean + lastP * valUB > targetM;
    };
    
    
    auto it = std::lower_bound(val, val + size, targetMean, f);
    
    
    double r = 0;
    if (true)
    {
      ing i = it - val;
      double lastP = 1.0 - psum[i - 1];
      double scaler = valUB / val[i];
      double partialMean = msum[i - 1] * scaler;
      r = (targetMean - lastP * valUB) / partialMean;
      r *= scaler;
    }
    
    
    return r;
    
  }
};







































