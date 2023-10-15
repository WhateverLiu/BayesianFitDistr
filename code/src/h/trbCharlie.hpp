
/*
struct BetaIntRaw
{
  bool earlyReturn;
  double a, b;
  double tgamma_a, tgamma_b;
  double r;
  double tgamma_a_plus_b;
  double apFinal, bpFinal;
  double lgamma_ap, lgamma_bp;
  
  
  bool notInt(double x)
  {
    return std::abs( x - int(x) ) > 1e-7 * std::max(1.0, std::abs(x) );
  };
  
  
  void reset(double a, double b)
  {
    this->a = a; this->b = b;
    tgamma_a = tgamma_b = 0; 
    if (b > 0)
    {
      tgamma_a = std::tgamma(a);
      tgamma_b = std::tgamma(b);
    }
    r = std::floor(-b);
    // earlyReturn = ! ( notInt(b) and a - r - 1 > 0 );
    earlyReturn = !notInt(b) or a - r - 1 <= 0 ;
    tgamma_a_plus_b = lgamma_ap = lgamma_bp = 0;
    if (earlyReturn) return;
    tgamma_a_plus_b = std::tgamma(a + b);
    apFinal = a - 1 - r; // apFinal = apFinal - 1 - r;
    bpFinal = b + 1 + r; // bpFinal = bpFinal + 1 + r;
    lgamma_ap = std::lgamma(apFinal);
    lgamma_bp = std::lgamma(bpFinal);
  }
  
  
  double operator()(double x, double x1m)
  {
    if (b > 0)
    {
      double Ix = (x > 0.5) ?
      R::pbeta(x1m, b, a, 0, 0) :
      R::pbeta(x,   a, b, 1, 0);
      return tgamma_a * tgamma_b * Ix;
    }
    
    
    // if ( ! ( notInt(b) and a - r - 1 > 0 ) ) { return 1e300; }
    if ( earlyReturn ) { return 1e300; }
    
    
    int i;
    double ap = a, bp = b;
    double lx = std::log(x);
    double lx1m = std::log(x1m);
    double x1 = std::exp(lx1m - lx);
    double c, tmp, sum, ratio;
    
    
    ap--;
    c = std::exp(ap * lx + bp * lx1m) / bp;
    sum = c;
    ratio = 1 / bp;
    bp++;
    
    
    for (i = 0; i < r; i++)
    {
      tmp = ap / bp;		
      c = tmp * (c * x1);	
      sum += c;
      ratio *= tmp;
      ap--;
      bp++;
    }
    
    
    double lIx = (x > 0.5) ? R::pbeta(x1m, bp, ap, 0, 1) : 
      R::pbeta(x, ap, bp, 1, 1);
    
    
    return -tgamma_a_plus_b * sum + (ratio * ap) * 
      std::exp(lgamma_ap + lgamma_bp + lIx);
  }
};
*/


// Regularized incomplete beta function.
// double rib(double x, double a, double b, int lowerTail = 1)
// {
//   return x <= 0.5 ? R::pbeta(x, a, b, lowerTail, 0) : 
//     R::pbeta( (0.5 - x) + 0.5 , b, a, 1 - lowerTail, 0);
// }


#include "bisecRoot.hpp"
#include "newtonRoot.hpp"


inline double softplus(double x) { 
  return std::log1p(std::exp(-std::abs(x))) + std::max(x, 0.0); 
}


struct TrB
{
  
  
  double a, b, c, d, lm1; // lm1 = E(min(X, 1)) where X ~ TrB
  double A, lnc, D, lnd; // lm1Gratio, lnlm1Gratio;
  // double tgamma_shape1_x_tgamma_shape3;
  // BetaIntRaw BIR;
  double lm1Gratio;
  
  
  double logPdf(double x) // After d is found.
  {
    x = std::max(x, 1e-300);
    double lnx = std::log(x);
    return A + (lnc - b * lnd) + (b - 1) * lnx + 
      D * softplus( c * (lnx - lnd) );
  }
  
  
  template <int RIBlib>
  double cdf(double q, int lowerTail = 1, int giveLog = 0)
  {
    if (q <= 0) return !lowerTail;
    double logvm, u;
    logvm = c * (lnd - std::log(q));
    u = std::exp(-softplus(logvm));
    if (u > 0.5)
    {
      double u1m = std::exp(-softplus(-logvm)); // log1pexp == softplus
      return ribetaf<RIBlib> (u1m, a / c, b / c, 1 - lowerTail, giveLog);
    }
    return ribetaf<RIBlib> (u, b / c, a / c, lowerTail, giveLog);
  }
  
  
  // Find quantile given a, b, c, d. d MUST have been solved before calling
  // this function. p should be in [0, 1)
  template <int RIBlib>
  double q(double p, int lowerTail = 0)
  {
    // We fixed it to Boost or R::pbeta since 
    // Numerical Recipes's solution is far from accurate.
    return d * std::pow(1.0 / ribetaf_inv<RIBlib % 2> (
        p, b/c, a/c, lowerTail, 0) - 1.0, -1.0 / c);
  }
  
  
  void resetTmpVars() // D, A, lnc, lm1Gratio.
  {
    D = -(a + b) / c;
    double lnGamma1 = std::lgamma(a / c);
    double lnGamma2 = std::lgamma(b / c);
    double lnGamma3 = lnGamma1 + lnGamma2;
    A = std::lgamma(-D) - lnGamma3;
    lnc = std::log(c);
    double lnlm1Gratio = std::lgamma( (b + 1) / c ) +
      std::lgamma( (a - 1) / c ) - lnGamma3;
    lm1Gratio = std::exp(lnlm1Gratio);
  } 
  
  
  void reset(double a, double b, double c, double d)
  {
    this->a = a; this->b = b; this->c = c; 
    resetTmpVars();
    this->d = d;
    lnd = std::log(d);
  }
  TrB(){}
  TrB(double a, double b, double c, double d){ reset(a, b, c, d); }
  
  
  template <int RIBlib>
  std::pair<double, double> computeL1m(
      double scale) // Given scale and fixed a, b, c.
  {
    double lnscale = std::log(scale);
    double clnscale = c * lnscale;
    double spclnscale = softplus(clnscale);
    double u = std::exp(-spclnscale);
    double br1, br2;
    double u1m = 0;
    if (u > 0.5)
    {
      u1m = std::exp(clnscale - spclnscale);
      br2 = ribetaf<RIBlib> (u1m, a / c, b / c, 1, 0);
    }
    else
    {
      br2 = ribetaf<RIBlib> (u, b / c, a / c, 0, 0);
    }
    if (br2 <= 5e-17) return std::pair<double, double> (
        scale * lm1Gratio, lm1Gratio);
    
    
    br1 = u > 0.5 ? ribetaf<RIBlib> (u1m, (a - 1) / c, (b + 1) / c, 0, 0) :
      ribetaf<RIBlib> (u, (b + 1) / c, (a - 1) / c, 1, 0);
    
    
    return std::pair<double, double> (
        scale * lm1Gratio * br1 + br2, lm1Gratio * br1);
    // Surprisingly, the derivative with respective to scale is just lm1Gratio * br1.
  }
  
  
  template <bool useNewton, int RIBlib>
  double reset(double a, double b, double c, 
               double lm1, double eps, int maxit)
  {
    this->a = a; this->b = b; this->c = c; this->lm1 = lm1;
    resetTmpVars();
    double initx = 1.0;
    if (!useNewton)
    {
      auto f = [&lm1, this](double scale)->double {
        return this->computeL1m<RIBlib> (scale).first - lm1; };
      d = bisecRoot(initx, f, eps, maxit);
    }
    else
    {
      auto f = [&lm1, this](double scale)->std::pair<double, double>
      {
        auto rst = this->computeL1m<RIBlib> (scale);
        rst.first -= lm1;
        return rst;
      };
      d = newtonRoot(initx, f, eps, maxit);
    }
    lnd = std::log(d);
    return d;
  }
  
  
  
  
  
  
  // ===========================================================================
  // Functions below are only relavant to discretization research, which in
  // end, is decided not to be used.
  // ===========================================================================
  std::vector<double> cnt;
  int supportSize;
  double *cdfGrad, infGrad;
  
  
  // xfirst > 0. The discretizer does not account for the probability at 0.
  template <bool useNewton, int RIBlib>
  void resetForDiscretization( // Tail probability will be loaded to the last point.
      double &xfirst, double &xlast, int size,
      double a, double b, double c, double lm1, 
      double eps, int maxit,
      double lastProbLowerBound,
      double *rstSupport,
      double *rstP
      ) // lastProbLowerBound equals 1e-10, for example.
  {
    supportSize = size;
    
    
    // Always use Newton to solve for initial d for the continous PDF.
    reset<true, RIBlib> (a, b, c, lm1, eps, maxit);
    double qval = q<RIBlib> (lastProbLowerBound, false);
    
    
    double xlastNew = std::min(xlast, qval);
    if (xlastNew != xlast)
    {
      xlast = xlastNew;
      xfirst = xlast / size;
    }
    if (useNewton) { cnt.resize(size); cdfGrad = &cnt[0]; }
    double dx = (xlast - xfirst) / (size - 1);
    double xmin = xfirst - dx / 2;
    for (int i = 0; i < size; ++i) rstSupport[i] = xmin + dx * i;
    if (useNewton) { infGrad = 0; }
  }
  
  
  
  
  template <int RIBlib>
  double giveTentativeMax(
      double a, double b, double c, double lm1, 
      double eps, int maxit, double lastProbLowerBound)
  {
    
    // Always use Newton to solve for initial d for the continous PDF.
    reset<true, RIBlib> (a, b, c, lm1, eps, maxit);
    return q<RIBlib> (lastProbLowerBound, false);
  }
  
  
  double diffCdf_d(double x)
  {
    x = std::max(x, 1e-300);
    double lnx = std::log(x);
    double lnx_lnd = lnx - lnd;
    double b_plus_c = b + c;
    double e = (D - 2) * softplus(c * lnx_lnd) + 
      b_plus_c * lnx_lnd + (2 * lnc - lnd - lnx) + A;
    return -std::exp(e);
  }
  
  
  // Given d, compute discretized mean and gradient.
  double psum;
  template <bool useNewton, int RIBlib, bool PfirstCoverFullDelta>
  std::pair<double, double> discretizedMean(
      double d, double *rstSupport, double *rstP, double p0)
  {
    std::pair<double, double> rst;
    this->d = d;
    lnd = std::log(d);
    for (int i = 0; i < supportSize; ++i)
    {
      rstP[i] = cdf<RIBlib> (rstSupport[i]);
      if (useNewton) cdfGrad[i] = diffCdf_d(rstSupport[i]);
    }
    
    
    double &meanval = rst.first, &meangrad = rst.second;
    meanval = 0; meangrad = 0; psum = 0;
    for (int i = 1; i < supportSize; ++i)
    {
      double mid = (rstSupport[i] + rstSupport[i - 1]) * 0.5;
      double pmass = 0;
      if (i == 1)
      {
        pmass = PfirstCoverFullDelta ? rstP[i] : rstP[i] - rstP[i - 1];
      }
      else pmass = rstP[i] - rstP[i - 1];
      meanval += mid * pmass;
      psum += pmass;
      if (useNewton) meangrad += mid * (cdfGrad[i] - cdfGrad[i - 1]);
    }
    
    
    if (true) // Deal with the last point.
    {
      int i = supportSize;
      double mid = (rstSupport[i - 1] - rstSupport[i - 2]) * 0.5 + rstSupport[i - 1];
      double pmass = 1.0 - rstP[i - 1];
      meanval += mid * pmass;
      psum += pmass;
      if (useNewton) meangrad += mid * (infGrad - cdfGrad[i - 1]);
    }
    
    
    meanval *= (1.0 - p0) / psum;
    if (useNewton) meangrad = (1.0 - p0) / psum;
    
    
    return rst;
  }
  
  
  
  
  // pmf[] is a given space of size doubles. Stores the result.
  template <bool useNewton, bool makeAndReturnSupport, int RIBlib,
            bool PfirstCoverFullDelta>
  void discretize(
      double &xfirst, double &xlast, int size,
      double a, double b, double c, double lm1, 
      double eps, int maxit,
      double lastProbLowerBound,
      double *support,
      double *pmf, double p0)
  {
    
    
    // Tail probability will be loaded to the last point.
    resetForDiscretization<useNewton, RIBlib> (
      xfirst, xlast, size, a, b, c, lm1,  eps, maxit, 
      lastProbLowerBound, support, pmf);
    
    
    auto objf = [this, support, pmf, p0](double scale)->std::pair<double, double>
    {
      auto rst = discretizedMean<useNewton, RIBlib, PfirstCoverFullDelta> (
        scale, support, pmf, p0);
      rst.first -= this->lm1;
      return rst;
    };
    
    
    if (useNewton) d = newtonRoot(d, objf, eps, maxit);
    else
    {
      auto objf2 = [&objf](double scale)->double { return objf(scale).first; };
      d = bisecRoot(d, objf2, eps, maxit);  
    }
    
    
    // Rcout << "d solved = " << d << "\n";
    
    
    lnd = std::log(d);
    
    
    // double invPsum = 1.0 / psum;
    double invPsum = (1.0 - p0) / psum;
    // At this moment, pmf still stores cdf.
    for (int i = 1; i < size; ++i)
    {
      if (i == 1)
      {
        pmf[i - 1] = PfirstCoverFullDelta ? pmf[i] * invPsum : 
          (pmf[i] - pmf[i - 1]) * invPsum;
      }
      else pmf[i - 1] = (pmf[i] - pmf[i - 1]) * invPsum;
    }

    
    pmf[size - 1] = (1.0 - pmf[size - 1]) * invPsum;
    
    
    if (makeAndReturnSupport)
    {
      double dx = support[1] - support[0], xmin = support[0] + dx / 2;
      for (int i = 0, iend = size - 1; i < iend; ++i) 
        support[i] = xmin + i * dx;
      support[size - 1] = xlast;  
    }
  }
  
  
  
  
  template <int RIBlib>
  double computeLimitedMean(double limit) // Provided that a, b, c, d have been computed. x is the limit.
  {
    // double scale = d;C
    // double lnscale = std::log(scale); // lnd exists.
    // double clnscale = c * lnscale;
    
    
    // Rcout << "a, b, c, d = ";
    // Rcout << a << " ";
    // Rcout << b << " ";
    // Rcout << c << " ";
    // Rcout << d << " ";
    
    
    double loglimit = std::log(limit);
    double clnscale = c * (lnd - loglimit);
    double spclnscale = softplus(clnscale);
    double u = std::exp(-spclnscale);
    double br1, br2;
    double u1m = 0;
    if (u > 0.5)
    {
      u1m = std::exp(clnscale - spclnscale);
      br2 = ribetaf<RIBlib> (u1m, a / c, b / c, 1, 0);
    }
    else
    {
      br2 = ribetaf<RIBlib> (u, b / c, a / c, 0, 0);
    }
    
    
    br2 *= limit;
    // if (br2 <= 5e-17) return std::pair<double, double> (
    //   scale * lm1Gratio, lm1Gratio);
    // if (br2 <= 5e-17) return d * lm1Gratio;
    
    
    br1 = u > 0.5 ? ribetaf<RIBlib> (u1m, (a - 1) / c, (b + 1) / c, 0, 0) :
      ribetaf<RIBlib> (u, (b + 1) / c, (a - 1) / c, 1, 0);
    
    
    // return std::pair<double, double> (
    //     scale * lm1Gratio * br1 + br2, lm1Gratio * br1);
    return d * lm1Gratio * br1 + br2;
    
  }
  
  
  
  
  // Given support, discretize the TrB onto the support
  // d has been solved at this moment, using the unbiased version.
  template <int RIBlib>
  void generatePMF(double *support, double *p, int size)
  {
#define lmtm computeLimitedMean<RIBlib>
    double *x = support;
    double h = x[1] - x[0];
    double psum = 0;
    p[0] = (lmtm(x[0]) - lmtm(x[1])) / h + 1.0 - cdf<RIBlib>(x[0]);
    psum += p[0];
    for (int i = 1, iend = size - 1; i < iend; ++i)
    { 
      p[i] = (2 * lmtm(x[i]) - lmtm(x[i - 1]) - 
        lmtm(x[i + 1])) / h;
      psum += p[i];
    }  
    p[size - 1] = (lmtm(x[size - 1]) - lmtm(x[size - 2])) / h - 
      1.0 + cdf<RIBlib>(x[size - 1]);
    psum += p[size - 1];
    double nr = 1.0 / psum;
    for (auto it = p, itend = p + size; it < itend; ++it) *it *= nr;
#undef lmtm 
  }
  
  
  // Given support, discretize the TrB onto the support
  // d has been solved at this moment.
  // Assume size >= 4, and the support is equispaced.
  template <int RIBlib>
  void unbiasedDiscretize(double *support, double *p, int size, double p0)
  {
    if (p0 >= 0)
    {
      support[0] = 0;
      p[0] = p0;
      support += 1;
      p += 1;
      size -= 1;
    }
    generatePMF<RIBlib>(support, p, size);
    if (p0 > 0)
    {
      double r = 1.0 - p0;
      for (auto pend = p + size; p < pend; ++p) *p *= r;
    }
  }
  
  
  template <int RIBlib>
  double unbiasedDiscretize( // Return d.
      double a, double b, double c, double lm1, 
      double eps, int maxit, double lastProbLowerBound,
      double *support, double *pmf, int size, double p0)
  {
    reset<true, RIBlib>(a, b, c, lm1, eps, maxit);
    if (p0 > 0) lastProbLowerBound /= 1.0 - p0;
    double maxval = std::min(1.0, q<RIBlib> (lastProbLowerBound));
    makeRegularGrid(0.0, maxval, support, size);
    unbiasedDiscretize<RIBlib> (support, pmf, size, p0);
    return d;
  }

  
  
  
};





















