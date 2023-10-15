// #include <Rcpp.h>
// using namespace Rcpp;


/*
template <bool useNewton, int RIBlib>
struct UnbiasedDiscretize
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
  
  
  void operator()(
      double *abc_lm1, double eps, int maxit, double tailProbLowerB,
      double *support, double *pmf, int size, double p0)
  {
    
    
    
    // unbiasedDiscretize(
    //   double a, double b, double c, double lm1, 
    //   double eps, int maxit, double lastProbLowerBound,
    //   double *support, double *pmf, int size, double p0)
    
    
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
*/


template <int RIBlib>
List makeDistrTableUnbiased(
    NumericMatrix paramMat, // NumericVector tentativeMax,
    NumericVector lastProbLowerBound, // size can be 1.
    NumericVector P0, // size can be 1.
    int size, // pmf size.
    double eps, 
    int maxit,  
    int maxCore)
{
  
  
  int Ndistr = paramMat.ncol();
  if (paramMat.nrow() != 4) stop("Input parameter has wrong shape.");
  if (P0.size() == 0) stop("Where is P0?");
  NumericMatrix Pmat(size + 2, Ndistr); // mean, max, probs.
  double *rst = &Pmat[0];
  double *abc_lm1 = &paramMat[0];
  int lastProbLowerBoundSize = lastProbLowerBound.size();
  int P0size = P0.size();
  CharlieThreadPool cp(maxCore);
  maxCore = cp.maxCore;
  int grainSize = Ndistr / (maxCore * maxCore * maxCore) + 1;
  std::vector<TrB> trbvec(maxCore);
  NumericVector dscale(Ndistr);
  std::vector<double> muErr(Ndistr);
  std::vector<std::vector<double> > supports(
      maxCore, std::vector<double> (size));
  
  
  cp.parFor(0, Ndistr, [&](std::size_t i, std::size_t t)->bool
  { 
    
    double *param = abc_lm1 + i * 4;
    double *entry = rst + i * (size + 2);
    entry[0] = param[3];
    dscale[i] = trbvec[t].unbiasedDiscretize<RIBlib> (
        param[0], param[1], param[2], param[3], eps, maxit, 
        lastProbLowerBound[i % lastProbLowerBoundSize],
        supports[t].data(), entry + 2, size, P0[i % P0size]);
    entry[1] = supports[t].back();
    
    
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




/*
template <int RIBlib>
List makeDistrTableUnbiased(
    NumericMatrix paramMat, // NumericVector tentativeMax,
    NumericVector lastProbLowerBound,
    NumericVector P0,
    int size, 
    double eps, 
    int maxit,  
    int maxCore)
{ 
  
  
  unBiasedDiscretize(
    paramMat,
    lastProbLowerBound,
    P0,
    size,
    eps, 
    maxit,  
    maxCore)
  
  
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
*/




// [[Rcpp::export]] 
List makeDistrTableUnbiased(
    NumericMatrix paramMat,
    NumericVector lastProbLowerBound,
    NumericVector P0,
    int size, 
    double eps, 
    int maxit,  
    int maxCore,
    String RIBlib = "Numerical Recipes")
{
  if (RIBlib == "Boost") return makeDistrTableUnbiased<0> (
    paramMat, lastProbLowerBound, P0, size, eps, maxit, maxCore);
  
  
  if (RIBlib == "R::pbeta") return makeDistrTableUnbiased<1> (
    paramMat, lastProbLowerBound, P0, size, eps, maxit, maxCore);
  
  
  if (RIBlib == "Numerical Recipes") return makeDistrTableUnbiased<2> (
    paramMat, lastProbLowerBound, P0, size, eps, maxit, maxCore);
  
  
  return List::create();
}




/*
template <int RIBlib>
List discretizeUnbiased(NumericVector abc_lm1, 
                double lastProbLowerBound,
                int size, 
                double eps, 
                int maxit)
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
*/



// [[Rcpp::export]]
List discretizeUnbiased(
    NumericMatrix abc_lm1, 
    NumericVector lastProbLowerBound,
    NumericVector P0,
    int size, 
    double eps, 
    int maxit,
    String RIBlib = "R::pbeta")
{
  
  List rst = makeDistrTableUnbiased(
    abc_lm1,
    lastProbLowerBound,
    P0,
    size, 
    eps, 
    maxit,  
    1,
    RIBlib);
  // List::create(Named("Ptable") = Pmat, Named("d") = dscale);
  
  
  NumericMatrix ps = rst["Ptable"];
  double mx = ps[1];
  NumericVector support(size);
  makeRegularGrid(0, mx, &support[0], size);
  return List::create(
    Named("val") = support, 
    Named("P") = NumericVector(ps.begin() + 2, ps.end()));
}














