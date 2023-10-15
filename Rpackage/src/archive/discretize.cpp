#include <Rcpp.h>
using namespace Rcpp;


template <bool useNewton, int RIBlib>
struct Discretize
{
  vec<double > val;
  TrB trb;
  void reset(int size)
  { 
    this->size = std::max(size, 3);
    val.resize(size);
  } 
  Discretize(){}
  Discretize(int size) { val.resize(size); }
  
  
  void operator()(double *abc_lm1,
                double tentMax, double tailProbLowerB,
                double *max, double *p, double *d,
                double eps, int maxit, double p0,
                bool PfirstCoverFullDelta)
  {
    double firstPoint = tentMax / val.size();
    
    
    if (PfirstCoverFullDelta) 
      trb.discretize<useNewton, false, RIBlib, true> (
          firstPoint, tentMax, val.size(), 
          abc_lm1[0], abc_lm1[1], abc_lm1[2], abc_lm1[3], eps, maxit, 
          tailProbLowerB, 
          &val[0], p, p0);
    else
      trb.discretize<useNewton, false, RIBlib, false> (
          firstPoint, tentMax, val.size(), 
          abc_lm1[0], abc_lm1[1], abc_lm1[2], abc_lm1[3], eps, maxit, 
          tailProbLowerB, 
          &val[0], p, p0);
    
    
    *d = trb.d;
    *max = val.back() + (val[1] - val[0]) / 2;
  }
};



template <bool useNewton, int RIBlib>
List discretize(
    // List discretizeNewton(
    NumericMatrix paramMat,
    NumericVector tentativeMax,
    NumericVector lastProbLowerBound,
    NumericVector P0,
    int size, 
    double eps, 
    int maxit,  
    int maxCore,
    bool PfirstCoverFullDelta)
{   
  int Ndistr = paramMat.ncol();
  if (paramMat.nrow() != 4) stop("Input parameter has wrong shape.");
  if (P0.size() == 0) stop("Where is P0?");
  NumericMatrix Pmat(size + 2, Ndistr); // mean, max, probs.
  double *rst = &Pmat[0];
  double *abc_lm1 = &paramMat[0];
  // bool p0notGiven = P0.size() == 0;
  bool singleMax = tentativeMax.size() == 1;
  bool singlePrLb = lastProbLowerBound.size() == 1;
  double *tm = &tentativeMax[0];
  double *sgp = &lastProbLowerBound[0];
  NumericVector dscale(Ndistr);
  double *dv = &dscale[0];
  // double zero = 0;
  // double *p0 = p0notGiven ? &zero : &P0[0];
  double *p0 = &P0[0];
  
  
  CharlieThreadPool cp(maxCore);
  maxCore = cp.maxCore;
  vec<Discretize<useNewton, RIBlib> > DZ(maxCore, size - 1);
  size = DZ[0].val.size() + 1;
  
  
  Discretize<useNewton, RIBlib> *dzv = &DZ[0];
  int grainSize = std::round(
    std::log(Ndistr) / std::log(std::max(2, maxCore))) + 1;
  
  
  vec<double> muErr(Ndistr);
  cp.parFor(0, Ndistr, [
              abc_lm1, rst, singleMax, singlePrLb, tm, sgp, 
              // p0notGiven, 
              p0, eps, maxit, size, dzv, dv, // &meanUnmatched, 
              &muErr , PfirstCoverFullDelta 
              // &whichThread
  ] 
  (std::size_t i, std::size_t t)->bool
  { 
    
    double *param = abc_lm1 + (i << 2);
    double *entry = rst + i * (size + 2);
    entry[0] = param[3];
    // double Pzero = p0notGiven ? 0 : p0[i];
    double Pzero = p0[i];
    // double Pbutzero = 1.0 - Pzero;
    // param[3] /= Pbutzero;
    double tmax = singleMax ? tm[0] : tm[i];
    double lb = singlePrLb ? sgp[0] : sgp[i];
    
    
    // lb /= Pbutzero;
    lb /= 1.0 - Pzero;
    
    
    dzv[t](param, tmax, lb, entry + 1, entry + 3, dv + i, eps, maxit, Pzero,
           PfirstCoverFullDelta);
    
    
    // param[3] = entry[0];
    
    
    // if (!p0notGiven)
    // {
    entry[2] = Pzero; 
    // double *pmf = entry + 3, *pmfend = entry + (size + 2);
    // for (; pmf < pmfend; ++pmf) *pmf *= Pbutzero;
    // }
    
    
    // Check if means are matched.
    if (true)
    { 
      double finalMean = 0;
      double *pmf = entry + 2;
      double dx = entry[1] / (size - 1);
      for (int u = 0; u < size; ++u) finalMean += pmf[u] * dx * u;
      muErr[i] = finalMean - entry[0];
    } 
    
    
    return false;
  }, grainSize,  
  [](std::size_t t)->bool{ return false; }, 
  [](std::size_t t)->bool{ return false; });
  
  
  for (int i = 0; i < Ndistr; ++i)
  { 
    if (std::abs(muErr[i]) > eps * 3)
      warning ( "MDR and mean differ by more than " +
        std::to_string(muErr[i]) );
  }
  
  
  return List::create(Named("Ptable") = Pmat, Named("d") = dscale);
} 




template <int RIBlib>
List makeDistrTable(NumericMatrix paramMat,
                    NumericVector tentativeMax,
                    NumericVector lastProbLowerBound,
                    NumericVector P0,
                    int size, 
                    double eps, 
                    int maxit,  
                    int maxCore,
                    bool useNewton = false,
                    bool PfirstCoverFullDelta = true)
{ 
  if (useNewton)
    return discretize<true, RIBlib> (
        paramMat,
        tentativeMax,
        lastProbLowerBound,
        P0,
        size, 
        eps, 
        maxit,  
        maxCore,
        PfirstCoverFullDelta);
   
  
  return discretize<false, RIBlib> (
      paramMat,
      tentativeMax,
      lastProbLowerBound,
      P0,
      size, 
      eps, 
      maxit,  
      maxCore,
      PfirstCoverFullDelta);
} 




// [[Rcpp::export]] 
List makeDistrTable(NumericMatrix paramMat,
                    NumericVector tentativeMax,
                    NumericVector lastProbLowerBound,
                    NumericVector P0,
                    int size, 
                    double eps, 
                    int maxit,  
                    int maxCore,
                    bool useNewton = false,
                    String RIBlib = "Numerical Recipes",
                    bool PfirstCoverFullDelta = true)
{
  if (RIBlib == "Boost") return makeDistrTable<0> (
    paramMat,tentativeMax, lastProbLowerBound, P0, size, eps, maxit, maxCore, useNewton,
    PfirstCoverFullDelta);
  
  
  if (RIBlib == "R::pbeta") return makeDistrTable<1> (
    paramMat,tentativeMax, lastProbLowerBound, P0, size, eps, maxit, maxCore, useNewton,
    PfirstCoverFullDelta);
  
  
  if (RIBlib == "Numerical Recipes") return makeDistrTable<2> (
    paramMat,tentativeMax, lastProbLowerBound, P0, size, eps, maxit, maxCore, useNewton,
    PfirstCoverFullDelta);
  
  
  return List::create();
}





template <int RIBlib>
List discretize(NumericVector abc_lm1, 
                double tentativeMax,
                double lastProbLowerBound,
                int size, 
                double eps, 
                int maxit,
                bool useNewton = false,
                bool PfirstCoverFullDelta = true)
{
  TrB trb;
  double firstPoint = tentativeMax / size;
  NumericVector val(size), p(size);
  if (useNewton)
  {   
    if (PfirstCoverFullDelta) 
      trb.discretize<true, true, RIBlib, true> (
          firstPoint, tentativeMax, size, 
          abc_lm1[0], abc_lm1[1], abc_lm1[2], abc_lm1[3], eps, maxit, 
          lastProbLowerBound, &val[0], &p[0], 0.0);
    else 
      trb.discretize<true, true, RIBlib, false> (
          firstPoint, tentativeMax, size, 
          abc_lm1[0], abc_lm1[1], abc_lm1[2], abc_lm1[3], eps, maxit, 
          lastProbLowerBound, &val[0], &p[0], 0.0); 
  }
  else
  {   
    if (PfirstCoverFullDelta) 
      trb.discretize<false, true, RIBlib, true> (
          firstPoint, tentativeMax, size, 
          abc_lm1[0], abc_lm1[1], abc_lm1[2], abc_lm1[3], eps, maxit, 
          lastProbLowerBound, &val[0], &p[0], 0.0);
    else
      trb.discretize<false, true, RIBlib, false> (
          firstPoint, tentativeMax, size, 
          abc_lm1[0], abc_lm1[1], abc_lm1[2], abc_lm1[3], eps, maxit, 
          lastProbLowerBound, &val[0], &p[0], 0.0);
  } 
  // Rcout << "trb.d = " << trb.d << "\n";
  
  
  NumericVector abcd(abc_lm1.begin(), abc_lm1.end());
  abcd[3] = trb.d;
  return List::create(
    Named("abcd") = abcd, 
    Named("distr") = DataFrame::create(Named("val") = val, Named("P") = p)
  );  
}




// [[Rcpp::export]]
List discretize(NumericVector abc_lm1, 
                double tentativeMax,
                double lastProbLowerBound,
                int size, 
                double eps, 
                int maxit,
                bool useNewton = false, 
                String RIBlib = "R::pbeta",
                bool PfirstCoverFullDelta = true)
{
  if (RIBlib == "Boost") return discretize<0> (
    abc_lm1, tentativeMax, lastProbLowerBound, size, eps, maxit, useNewton,
    PfirstCoverFullDelta);
  
  
  if (RIBlib == "R::pbeta") return discretize<1> (
    abc_lm1, tentativeMax, lastProbLowerBound, size, eps, maxit, useNewton,
    PfirstCoverFullDelta);
  
  
  if (RIBlib == "Numerical Recipes") return discretize<2> (
    abc_lm1, tentativeMax, lastProbLowerBound, size, eps, maxit, useNewton,
    PfirstCoverFullDelta);
  
  
  stop("Method of computing Regularized Incomplete Beta function is not implemented.");
  
  
  return List::create();
}














