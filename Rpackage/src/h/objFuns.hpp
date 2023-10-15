#pragma once
#include "trbCharlie.hpp"
#include "distance.hpp"


// distanceFun == "llh", 0: use the likelihood criterion
// distanceFun == "dlh", 1: use the likelihood criterion over discrete distribution.
// distanceFun == "kmg", 2: use the Komogorov distance criterion.
// distanceFun == "euc", 3: use the Euclidean distance criterion


template <int method, int RIBlib>
struct ObjFun
{
  int maxit, size;
  double lm1, eps;
  double *xst, *pmf;
  TrB trb;
  
  
  Negllh D0;  
  NegllhDiscrete D1;
  Kolmogorov D2;
  EucCDF D3; 
  UnbiasedNegllhDiscrete D4;
  UnbiasedKolmogorov D5;
  UnbiasedEucCDF D6;
  
  
  // FitObj(){ gradEps = 1.49011611938477e-8; }
  
  
  void reset(double *xst, double *pmf, int size,
             double lm1, double eps, int maxit)
  { 
    this->xst = xst;
    this->pmf = pmf;
    this->size = size;
    this->lm1 = lm1;
    this->eps = eps;
    this->maxit = maxit;
  } 
  
  
  double f(double *abc)
  { 
    
    // Compute d.
    trb.reset<true, RIBlib> (abc[0], abc[1], abc[2], lm1, eps, maxit);
    
    
    auto trbPtr = &trb;
    auto logPdf = [trbPtr](double t)->double { return trbPtr->logPdf(t); };
    auto cdf    = [trbPtr](double q)->double { return trbPtr->cdf<RIBlib> (q); };
    auto PMFgen = [trbPtr](double *x, double *p, int size)->void { 
      trbPtr->generatePMF<RIBlib> (x, p, size); };
    
    
    double rst = 0;
    if (method == 0)
    { 
      double logTailProb = trb.cdf<RIBlib> (1.0, 0, 1);
      rst = D0(xst, pmf, size, logPdf, logTailProb);
    } 
    else if (method == 1) rst = D1(xst, pmf, size, cdf);
    else if (method == 2) rst = D2(xst, pmf, size, cdf);
    else if (method == 3) rst = D3(xst, pmf, size, cdf);
    else if (method == 4) rst = D4(xst, pmf, size, PMFgen);
    else if (method == 5) rst = D5(xst, pmf, size, PMFgen);
    else if (method == 6) rst = D6(xst, pmf, size, PMFgen);
    return rst;
  } 
  
  
  double operator()(double *abc, double *&grad, int dim = 3)
  { 
    grad = nullptr;
    return f(abc);
  } 
};


