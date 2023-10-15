

// 1: Exponential scaling, 2: Gaussian scaling.
template <int method> // 
struct ExpScaleProb
{
  // x in (0, 1]
  std::vector<double> cntr;
  void operator()(
      double *x, double *p, int size, double targetMean, double *rstP)
  {
    cntr.resize(size);
    double *lx = cntr.data();
    if (method == 1)
    {
      for (int i = 0; i < size; ++i) lx[i] = std::log(x[i]);  
    }
    else
    {
      for (int i = 0; i < size; ++i) 
      {
        lx[i] = std::log(x[i]);   
        lx[i] *= -lx[i];
      }
    }
    
    
    auto f = [&](double u)->double 
    {
      
    };
    
  }
  
  
};


