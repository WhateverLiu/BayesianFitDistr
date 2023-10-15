
// #############################################################################
// Distance between a PMF and a continous distribution. The continuous
//   distribution is given by either its CDF or log PDF.
// #############################################################################


// Less Than But Close To 1
inline bool ltbct1(double x) { return x > 1 - 1e-10; }


// =============================================================================
// Essentially cross entropy. Mean log-likelihood.
// =============================================================================
struct Negllh
{
  template <typename logPDF>
  double operator()(
      double *x, double *pmf, int size, 
      logPDF &lpdf, double logTailProb) 
  {
    auto &p = pmf;
    double S = 0;
    int i = 0;
    for (int iend = size - 1; i < iend; ++i) 
      S += p[i] == 0 ? 0 : p[i] * lpdf(x[i]);
    if ( ltbct1(x[i]) ) S += p[i] * logTailProb;
    else S += p[i] * lpdf(x[i]);
    double rst = -S / size;
    return rst;
  }
  
  
  Negllh(){}
  double logTailProb;
  Negllh(double logTailProb): logTailProb(logTailProb) {}
  template <typename logPDF>
  double operator()(
      double *x, double *pmf, int size, logPDF &lpdf)
  {
    return this->operator()<logPDF>(x, pmf, size, lpdf, logTailProb);
  }
    
  
};


/* Unoptimized implementation
 // =============================================================================
 // Essentially cross entropy. Mean log-likelihood.
 // =============================================================================
 struct NegllhDiscrete
 {
 template <typename CDF> // x and pmf are equispaced. Assume x in [0, 1]
 double operator()(double *x, double *pmf, int size, CDF &cdf)
 {
 double gap = x[1] - x[0];
 double xprior = x[0] - gap / 2;
 double xpriorCdf = cdf(xprior);
 truePmf.resize(size);
 double psum = 0;
 int i = 0;
 double xnow, xnowCdf, pmass;
 xnow = xnowCdf = pmass = 0;
 for (int iend = size - 1; i < iend; ++i)
 {
 xnow = xprior + gap;
 xnowCdf = cdf(xnow);
 pmass = xnowCdf - xpriorCdf;
 truePmf[i] = pmass;
 psum += pmass;
 xprior = xnow;
 xpriorCdf = xnowCdf;
 }
 if ( x[i] - 1 > -1e-10 ) // If the last support point is very close to 1.
 {
 pmass = 1.0 - cdf(xprior);
 }
 else
 {
 xnow = xprior + gap;
 xnowCdf = cdf(xnow);
 pmass = xnowCdf - xpriorCdf;
 }
 truePmf[i] = pmass;
 psum += pmass;
 for (i = 0; i < size; ++i) Rcout << truePmf[i] << ", ";
 Rcout << "\n";
 double r = 1.0 / psum;
 for (i = 0; i < size; ++i) truePmf[i] *= r;
 double S = 0;
 for (i = 0; i < size; ++i) S += std::log(truePmf[i]) * pmf[i];
 return -S / size;
 }
 
 
 vec<double> truePmf;
 };
 */


// =============================================================================
// Essentially cross entropy. Mean log-likelihood.
// Since the empirical distribution does not change, this is also equivalent to
// Kullback.
// =============================================================================
struct NegllhDiscrete
{
  template <typename CDF> // x and pmf are equispaced. Assume x in [0, 1]
  double operator()(double *x, double *pmf, int size, CDF &cdf)
  {
    double gap = x[1] - x[0];
    double xprior = x[0] - gap / 2;
    double xpriorCdf = cdf(xprior);
    double psum = 0;
    double xnow, xnowCdf, pmass;
    xnow = xnowCdf = pmass = 0;
    double S = 0;
    int i = 0;
    for (int iend = size - 1; i < iend; ++i)
    {
      xnow = xprior + gap;
      xnowCdf = cdf(xnow);
      pmass = xnowCdf - xpriorCdf;
      S += std::log(pmass) * pmf[i];
      psum += pmass;
      xprior = xnow;
      xpriorCdf = xnowCdf;
    }
    if ( ltbct1(x[i]) ) // If the last support point is very close to 1.
    {
      pmass = 1.0 - cdf(xprior);
    }
    else
    {
      xnow = xprior + gap;
      xnowCdf = cdf(xnow);
      pmass = xnowCdf - xpriorCdf;
    }
    S += std::log(pmass) * pmf[i];
    psum += pmass;
    double r = 1.0 / psum;
    S += std::log(r);
    return -S / size;
  }
};




// =============================================================================
// Max distance between CDFs. Here the max operator is Boltzmann Smooth Max so
//   the function is differentiable.
// =============================================================================

struct KolmogorovUnadjusted 
{ 
  std::vector<double> d;
  template <typename CDF>
  double operator()(double *x, double *pmf, int size, CDF &cdf)
  {
    auto &p = pmf;
    d.resize(size);
    double psum = 0, dmax = 0;
    int i = 0;
    for (int iend = size - 1; i < iend; ++i)
    {
      psum += p[i];
      d[i] = std::abs(psum - cdf(x[i]));
      dmax = std::max(dmax, d[i]);
    }
    if ( ltbct1(x[i]) )
    {
      psum += p[i];
      d[i] = std::abs(psum - 1.0  );
      dmax = std::max(dmax, d[i]);
    }
    else
    {
      psum += p[i];
      d[i] = std::abs(psum - cdf(x[i]));
      dmax = std::max(dmax, d[i]);
    }
    if (dmax <= 0) return 0;
    constexpr const double alpha = 16;
    double numtr = 0, divider = 0;
    double dmaxInv = 1.0 / dmax;
    for (int i = 0; i < size; ++i)
    {
      double x = d[i] * dmaxInv;
      double e = std::exp(alpha * x);
      divider += e;
      numtr += x * e;
    }
    return numtr / divider * dmax;
  }
  
  
};




struct Kolmogorov 
{ 
  std::vector<double> d;
  template <typename CDF>
  double operator()(double *x, double *pmf, int size, CDF &cdf)
  { 
    auto &p = pmf;
    d.resize(size);
    double psum = 0, dmax = 0, delta = x[1] - x[0];
    for (int i = size - 1; i >= 0; --i)
    { 
      psum += p[i];
      d[i] = std::abs( 1.0 - cdf(x[i] - delta / 2) - psum );
      dmax = std::max(dmax, d[i]);
    }
    
    
    if (dmax <= 0) return 0;
    constexpr const double alpha = 16;
    double numtr = 0, divider = 0;
    double dmaxInv = 1.0 / dmax;
    for (int i = 0; i < size; ++i)
    { 
      double x = d[i] * dmaxInv;
      double e = std::exp(alpha * x);
      divider += e;
      numtr += x * e;
    } 
    return numtr / divider * dmax;
  } 
  
  
}; 






// =============================================================================
// Mean Euclidean distance between CDFs.
// =============================================================================
/*
struct EucCDF
{ 
  template <typename CDF>
  double operator()(double *x, double *pmf, int size, CDF &cdf)
  {
    const auto &p = pmf;
    double psum = 0;
    double rst = 0;
    int i = 0;
    for (int iend = size - 1; i < iend; ++i)
    { 
      psum += p[i];
      double d = psum - cdf(x[i]);
      rst += d * d;
    } 
    if ( ltbct1(x[i]) )
    { 
      psum += p[i];
      double d = psum - 1.0;
      rst += d * d;
    } 
    else
    { 
      psum += p[i];
      double d = psum - cdf(x[i]);
      rst += d * d;
    } 
    return std::sqrt(rst / size);
  }
};  
*/




struct EucCDF
{ 
  template <typename CDF>
  double operator()(double *x, double *pmf, int size, CDF &cdf)
  {
    auto &p = pmf;
    double psum = 0;
    double delta = x[1] - x[0];
    double dsum = 0;
    for (int i = size - 1; i >= 0; --i)
    { 
      psum += p[i];
      double u = 1.0 - cdf(x[i] - delta / 2) - psum;
      dsum += u * u;
    } 
    return std::sqrt(dsum / size);
    
    
    /*
    if (dmax <= 0) return 0;
    const auto &p = pmf;
    double psum = 0;
    double rst = 0;
    int i = 0;
    for (int iend = size - 1; i < iend; ++i)
    {  
      psum += p[i];
      double d = psum - cdf(x[i]);
      rst += d * d;
    }  
    if ( ltbct1(x[i]) )
    {  
      psum += p[i];
      double d = psum - 1.0;
      rst += d * d;
    }  
    else
    {  
      psum += p[i];
      double d = psum - cdf(x[i]);
      rst += d * d;
    }  
    return std::sqrt(rst / size);
    */
  } 
};  



struct UnbiasedNegllhDiscrete
{
  template <typename PMFgenerator> // x and pmf are equispaced. Assume x in [0, 1]
  double operator()(double *x, double *pmf, int size, PMFgenerator &Pgen)
  {
    
    // auto &p = pmf;
    d.resize(size);
    Pgen(x, d.data(), size);
    
    
    double nr = 1.0/ std::accumulate(d.begin(), d.end(), 0.0);
    for (auto it = d.begin(); it < d.end(); ++it) *it *= nr;
    
    
    double *pp = d.data(); // Generate the other PMF.
    double S = 0;
    for (int i = 0; i < size; ++i)
      S += std::log(pp[i]) * pmf[i];
    return -S / size;
  }
  vec<double> d;
};




struct UnbiasedKolmogorov
{
  // Unbiased distance between survival functions
  template <typename PMFgenerator>
  double operator()(double *x, double *pmf, int size, PMFgenerator &Pgen)
  {
    // auto &p = pmf;
    d.resize(size);
    Pgen(x, d.data(), size);
    double *pp = d.data(); // Generate the other PMF.
    
    
    double nr = 1.0/ std::accumulate(d.begin(), d.end(), 0.0);
    for (auto it = d.begin(); it < d.end(); ++it) *it *= nr;
    
    
    double pmfSum = pmf[size - 1];
    double ppSum = pp[size - 1];
    double dmax = std::abs(ppSum - pmfSum);
    for (int i = size - 2; i >= 0; --i)
    {
      pmfSum += pmf[i];
      ppSum += pp[i];
      double diff = std::abs(ppSum - pmfSum);
      dmax = std::max(diff, dmax);
      d[i] = diff;
    }
    
    
    if (dmax <= 0) return 0;
    constexpr const double alpha = 16;
    double numtr = 0, divider = 0;
    double dmaxInv = 1.0 / dmax;
    for (int i = 0; i < size; ++i)
    { 
      double x = d[i] * dmaxInv;
      double e = std::exp(alpha * x);
      divider += e;
      numtr += x * e;
    } 
    return numtr / divider * dmax;
  } 
  
  
  vec<double> d;
};




struct UnbiasedEucCDF
{
  // Unbiased distance between survival functions
  template <typename PMFgenerator>
  double operator()(double *x, double *pmf, int size, PMFgenerator &Pgen)
  { 
    // auto &p = pmf;
    d.resize(size);
    Pgen(x, d.data(), size);
    
    
    double nr = 1.0/ std::accumulate(d.begin(), d.end(), 0.0);
    for (auto it = d.begin(); it < d.end(); ++it) *it *= nr;
    
    
    double *pp = d.data(); // Generate the other PMF.
    
    
    double pmfSum = pmf[size - 1];
    double ppSum = pp[size - 1];
    double diff = ppSum - pmfSum;
    double dsum = diff * diff;
    for (int i = size - 2; i >= 0; --i)
    { 
      pmfSum += pmf[i];
      ppSum += pp[i];
      diff = ppSum - pmfSum;
      dsum += diff * diff;
    } 
    return std::sqrt(dsum / size);
  }  
  
  
  vec<double> d;
};



