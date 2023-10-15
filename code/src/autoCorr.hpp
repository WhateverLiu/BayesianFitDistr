

// [[Rcpp::export]]
DataFrame autoCorr(NumericVector x, int maxCore = 1000)
{
  const int xsize = x.size();
  if (xsize <= 3) stop("Data size too small.");
  
  
  struct mm { double s, ss; }; // sum and sum of squares.
  std::vector<mm> C(xsize + 1);
  C[0].s = C[0].ss = 0;
  mm *cs = C.data() + 1;
  for (int i = 0; i < xsize; ++i)
  {
    cs[i].s = cs[i - 1].s + x[i];
    cs[i].ss = cs[i - 1].ss + x[i] * x[i];
  }
  
  
  auto rangeSum = [&](int first, int last)->double 
    { return cs[last].s - cs[first - 1].s; };
  auto rangeSS = [&](int first, int last)->double 
    { return cs[last].ss - cs[first - 1].ss; };
  
    
  double *x_ = &x[0];
  auto acorlagk = [&](int k)->double
  {
    double inv_size = 1.0 / (xsize - k);
    
    
    double EA = rangeSum(0, xsize - k - 1) * inv_size;
    double EA2 = rangeSS(0, xsize - k - 1) * inv_size;
    double varA = EA2 - EA * EA;
    
    
    double EB = rangeSum(k, xsize - 1) * inv_size;
    double EB2 = rangeSS(k, xsize - 1) * inv_size;
    double varB = EB2 - EB * EB;
    
    
    double EAB = std::inner_product(x_, x_ + xsize - k, x_ + k, 0.0) * inv_size;
    return (EAB - EA * EB) / std::sqrt(varA * varB);
  };
  
  
  CharlieThreadPool cp(maxCore);
  maxCore = cp.maxCore;
  int Njob = xsize / 2 + 1 - 1;
  int grainSize = Njob / (maxCore * maxCore * maxCore) + 1;
  NumericVector corr(Njob);
  double *corr_ = &corr[0];
  
  
  // i is the lag.
  cp.parFor(1, xsize / 2 + 1, [&](std::size_t i, std::size_t t)->bool
  {
    corr_[i - 1] = acorlagk(i);
    return false;
  }, grainSize, 
  [](std::size_t t)->bool { return false; },
  [](std::size_t t)->bool { return false; });
  
  
  IntegerVector lag(Njob);
  std::iota(lag.begin(), lag.end(), 1);
  return DataFrame::create(Named("lag") = lag, Named("corr") = corr);
}



















