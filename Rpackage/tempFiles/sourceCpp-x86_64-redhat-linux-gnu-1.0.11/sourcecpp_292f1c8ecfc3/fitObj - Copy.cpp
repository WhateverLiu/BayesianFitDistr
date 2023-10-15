// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "h/objFuns.hpp"
#include "lbfgs/LBFGSBcharlie.hpp"
#include "h/charlieThreadPool2.hpp"




// OPT is the type of the optimizer.
template <int method, int RIBlib>
std::pair<double, int> solveShapeParamsGivenMethodLBFGS(
    double *abc_lm1, double *abcLB, double *abcUB,
    double *val, double *P, int size, // Empirical PMF.
    double scaleEps, int scaleMaxit,
    // =========================================================================
    // Optimizer options.
    // =========================================================================
    LBFGSBcharlie &opt,
    int max_iterations,
    double hgrad,
    bool centralDiff,
    double epsilon,
    double epsilon_rel,
    int m,
    int past,
    double delta,
    int max_submin,
    int max_linesearch,
    double min_step,
    double max_step,
    double ftol,
    double wolfe
  // ===========================================================================
) 
{ 
  
  
  ObjFun<method, RIBlib> objf;
  objf.reset(val, P, size, abc_lm1[3], scaleEps, scaleMaxit);
  
   
  auto rst = opt (
    abc_lm1, 3, objf, 
    max_iterations, 
    abcLB, 
    abcUB,
    hgrad,
    centralDiff,
    epsilon,
    epsilon_rel,
    m,
    past,
    delta,
    max_submin,
    max_linesearch,
    min_step,
    max_step,
    ftol,
    wolfe
  );
  
  
  return rst;
} 




// distanceFun == "likelihood",         0: use the likelihood criterion
// distanceFun == "likelihoodDiscrete", 1: use the likelihood criterion over discrete distribution.
// distanceFun == "Komogorov",          2: use the Komogorov distance criterion.
// distanceFun == "EuclideanCDF",       3: use the Euclidean distance criterion




template <int RIBlib>
List LBFGSBtrbFit(
    NumericVector abc_lm1,
    DataFrame empDistr,
    NumericVector abcLB,
    NumericVector abcUB,
    double scaleEps,
    int scaleMaxit,
    String distanceFun,
    int max_iterations,
    double hgrad,
    bool centralDiff,
    int m,
    double epsilon,
    double epsilon_rel,
    int past,
    double delta,
    int max_submin,
    int max_linesearch,
    double min_step,
    double max_step,
    double ftol,
    double wolfe
  )
{
  
  
  LBFGSBcharlie opt;
  NumericVector init(abc_lm1.begin(), abc_lm1.end());
  double *abc_lm1_ = &init[0], *abcLB_ = nullptr, *abcUB_ = nullptr;
  if (abcLB.size() != 0) abcLB_ = &abcLB[0];
  if (abcUB.size() != 0) abcUB_ = &abcUB[0];
  NumericVector val = empDistr[0], P = empDistr[1];
  double *val_ = &val[0], *P_ = &P[0];
  
  
  std::pair<double, int> rst(1e300, 0);
  
  
  if (false)
  {
    
    
    // double abc_lm1[4] = { abc_[0], abc_[1], abc_[2], lm1_[0] };
    for (int u = 0; u < 4; ++u) Rcout << init[u] << ", ";
    Rcout << "\n";
    
    
    for (int u = 0; u < 3; ++u) Rcout << abcLB_[u] << ", ";
    Rcout << "\n";
    
    
    for (int u = 0; u < 3; ++u) Rcout << abcUB_[u] << ", ";
    Rcout << "\n";
    
    
    Rcout << "method = " << 0 << ", ";
    Rcout << "RIBlib = " << RIBlib << "\n";
    
    
    for (int u = 0; u < val.size(); ++u)
      Rcout << val_[u] << ", ";
    Rcout << "\n";
    for (int u = 0; u < val.size(); ++u)
      Rcout << P_[u] << ", ";
    Rcout << "\n\n";
    
    
    Rcout << ", " << scaleEps << ", " << scaleMaxit << ", " << max_iterations << 
      ", " << hgrad << ", " << centralDiff << ", " << epsilon << ", " << 
        epsilon_rel << ", " << m << ", " << past << ", " << delta << ", " << 
          max_submin << ", " << max_linesearch << ", " << min_step << ", " << 
            max_step << ", " << ftol << ", " << wolfe << ", " << "\n";
  }
  
  


  if (distanceFun == "likelihood") rst = solveShapeParamsGivenMethodLBFGS<0, RIBlib> (
    abc_lm1_, abcLB_, abcUB_, val_, P_, val.size(), scaleEps,
    scaleMaxit, opt, max_iterations, hgrad, centralDiff, epsilon, 
    epsilon_rel, m, past, delta, max_submin, max_linesearch, min_step,
    max_step, ftol, wolfe
  );


  else if (distanceFun == "likelihoodDiscrete") 
    rst =  solveShapeParamsGivenMethodLBFGS<1, RIBlib> (
        abc_lm1_, abcLB_, abcUB_, val_, P_, val.size(), scaleEps,
        scaleMaxit, opt, max_iterations, hgrad, centralDiff, epsilon, 
        epsilon_rel, m, past, delta, max_submin, max_linesearch, min_step,
        max_step, ftol, wolfe
    );
  
  
  else if (distanceFun == "Komogorov") 
    rst =  solveShapeParamsGivenMethodLBFGS<2, RIBlib> (
        abc_lm1_, abcLB_, abcUB_, val_, P_, val.size(), scaleEps,
        scaleMaxit, opt, max_iterations, hgrad, centralDiff, epsilon, 
        epsilon_rel, m, past, delta, max_submin, max_linesearch, min_step,
        max_step, ftol, wolfe
    );
  
  
  else if (distanceFun == "EuclideanCDF") 
    rst =  solveShapeParamsGivenMethodLBFGS<3, RIBlib> (
        abc_lm1_, abcLB_, abcUB_, val_, P_, val.size(), scaleEps,
        scaleMaxit, opt, max_iterations, hgrad, centralDiff, epsilon, 
        epsilon_rel, m, past, delta, max_submin, max_linesearch, min_step,
        max_step, ftol, wolfe
    );
  
  
  else if (distanceFun == "unbiasedLikelihoodDiscrete") 
    rst =  solveShapeParamsGivenMethodLBFGS<4, RIBlib> (
      abc_lm1_, abcLB_, abcUB_, val_, P_, val.size(), scaleEps,
      scaleMaxit, opt, max_iterations, hgrad, centralDiff, epsilon, 
      epsilon_rel, m, past, delta, max_submin, max_linesearch, min_step,
      max_step, ftol, wolfe
    );
  
  
  else if (distanceFun == "unbiasedKomogorov") 
    rst =  solveShapeParamsGivenMethodLBFGS<5, RIBlib> (
      abc_lm1_, abcLB_, abcUB_, val_, P_, val.size(), scaleEps,
      scaleMaxit, opt, max_iterations, hgrad, centralDiff, epsilon, 
      epsilon_rel, m, past, delta, max_submin, max_linesearch, min_step,
      max_step, ftol, wolfe
    );
  
  
  else if (distanceFun == "unbiasedEuclideanCDF") 
    rst =  solveShapeParamsGivenMethodLBFGS<6, RIBlib> (
      abc_lm1_, abcLB_, abcUB_, val_, P_, val.size(), scaleEps,
      scaleMaxit, opt, max_iterations, hgrad, centralDiff, epsilon, 
      epsilon_rel, m, past, delta, max_submin, max_linesearch, min_step,
      max_step, ftol, wolfe
    );
  
  
  else stop("method not implemented");
  
  
  return List::create(Named("param") = init, Named("fval") = rst.first, 
               Named("niter") = rst.second);
}




// [[Rcpp::export]]
List LBFGSBtrbFit(
    NumericVector abc_lm1,
    DataFrame empDistr,
    NumericVector abcLB = NumericVector(0),
    NumericVector abcUB = NumericVector(0),
    double scaleEps = 1e-8, 
    int scaleMaxit = 100,
    String distanceFun = "likelihood",
    int max_iterations = 100,
    String RIBlib = "R::pbeta",
    double hgrad = 0,
    bool centralDiff = true,
    int m = 6,
    double epsilon = 1e-5,
    double epsilon_rel = 1e-5,
    int past = 1,
    double delta = 1e-10,
    int max_submin = 10,
    int max_linesearch = 20,
    double min_step = 1e-20,
    double max_step = 1e+20,
    double ftol = 1e-4,
    double wolfe = 0.9
)
{
  
  if (RIBlib == "Boost") return LBFGSBtrbFit<0>(
    abc_lm1, empDistr, abcLB, abcUB, scaleEps, scaleMaxit, distanceFun, max_iterations,
    hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta, max_submin,
    max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  if (RIBlib == "R::pbeta") return LBFGSBtrbFit<1>(
    abc_lm1, empDistr, abcLB, abcUB, scaleEps, scaleMaxit, distanceFun, max_iterations,
    hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta, max_submin,
    max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  if (RIBlib == "Numerical Recipes") return LBFGSBtrbFit<2>(
    abc_lm1, empDistr, abcLB, abcUB, scaleEps, scaleMaxit, distanceFun, max_iterations,
    hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta, max_submin,
    max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  stop("Method of computing Regularized Incomplete Beta function is not implemented.");
  
  
  return List::create();
}
  





// abc.ncol() == the number of initial parameters to test.
// abc.nrow() == 3.
// For every PMF in empDistrList, we test abc.ncol() sets of different parameters.
// lm1 == the limited means for all empirical PMFs.
// lm1.size() should be equal to empDistrList.size().
template <int method, int RIBlib>
List LBFGSBtrbFitListGivenMethodPlain(
    NumericMatrix abc,
    NumericVector lm1,
    List empDistrList,
    NumericVector abcLB,
    NumericVector abcUB,
    double scaleEps,
    int scaleMaxit,
    int max_iterations,
    int maxCore,
    double hgrad,
    bool centralDiff,
    int m,
    double epsilon,
    double epsilon_rel,
    int past,
    double delta,
    int max_submin,
    int max_linesearch,
    double min_step,
    double max_step,
    double ftol,
    double wolfe
)
{
  
  if (abc.nrow() != 3) stop("abc's number of rows != 3");
  
  
  if (lm1.size() != empDistrList.size())
    stop("empDistrList's size unequal to lm1's size.");
  double *lm1_ = &lm1[0];
  
  
  int nparam = abc.ncol();
  double *abc_ = &abc[0];
  double *abcLB_ = nullptr, *abcUB_ = nullptr;
  if (abcLB.size() != 0) abcLB_ = &abcLB[0];
  if (abcUB.size() != 0) abcUB_ = &abcUB[0];
  
  
  struct DistrPtr { int size; double *val, *P; };
  std::vector<DistrPtr> distrPtrs(empDistrList.size());
  for (int i = 0, iend = distrPtrs.size(); i < iend; ++i)
  {
    List distr = empDistrList[i];
    NumericVector val = distr[0], P = distr[1];
    distrPtrs[i].size = val.size();
    distrPtrs[i].val = &val[0];
    distrPtrs[i].P = &P[0];
  }
  DistrPtr *distrPtrs_ = &distrPtrs[0];
  
  
  NumericMatrix rstParam(4, distrPtrs.size());
  double *rstParam_ = &rstParam[0];
  
  
  NumericVector fval(distrPtrs.size());
  double *fval_ = &fval[0];
  
  
  IntegerVector niter(fval.size());
  int *niter_ = &niter[0];
  
  
  maxCore = std::min(nparam, maxCore);
  CharlieThreadPool cp(std::move(maxCore));
  // maxCore = cp.maxCore;
  std::vector<LBFGSBcharlie> opt(maxCore);
  LBFGSBcharlie *opt_ = &opt[0];
  
  
  // int grainSize = (fval.size() + 0.0) / (maxCore * maxCore * maxCore) + 1;
  int grainSize = fval.size() / (maxCore * maxCore * maxCore) + 1;
  // Rcout << grainSize << "\n";
  
  
  if (false)
  {
    double abc_lm1[4] = { abc_[0], abc_[1], abc_[2], lm1_[0] };
    for (int u = 0; u < 4; ++u) Rcout << abc_lm1[u] << ", ";
    Rcout << "\n";
    
    
    for (int u = 0; u < 3; ++u) Rcout << abcLB_[u] << ", ";
    Rcout << "\n";
    
    
    for (int u = 0; u < 3; ++u) Rcout << abcUB_[u] << ", ";
    Rcout << "\n";
    
    
    Rcout << "method = " << method << ", ";
    Rcout << "RIBlib = " << RIBlib << "\n";
    
    
    for (int u = 0; u < distrPtrs_[0].size; ++u)
      Rcout << distrPtrs_[0].val[u] << ", ";
    Rcout << "\n";
    for (int u = 0; u < distrPtrs_[0].size; ++u)
      Rcout << distrPtrs_[0].P[u] << ", ";
    Rcout << "\n";
    
    
    Rcout << ", " << scaleEps << ", " << scaleMaxit << ", " << max_iterations << ", " << hgrad << ", " << centralDiff << ", " << 
      epsilon << ", " << epsilon_rel << ", " << m << ", " << past << ", " << delta << ", " << max_submin << ", " << max_linesearch << ", " << 
        min_step << ", " << max_step << ", " << ftol << ", " << wolfe << ", " << "\n";
    
    
    auto rst = solveShapeParamsGivenMethodLBFGS<method, RIBlib> (
      abc_lm1, abcLB_, abcUB_, 
      distrPtrs_[0].val, distrPtrs_[0].P, distrPtrs_[0].size, 
      scaleEps, scaleMaxit, opt_[0], max_iterations, hgrad, centralDiff, 
      epsilon, epsilon_rel, m, past, delta, max_submin, max_linesearch, 
      min_step, max_step, ftol, wolfe
    );
    
    
    Rcout << "rst = " << rst.first << ", " << rst.second << "\n";
  }
  
  
  auto runner = [&](std::size_t i, std::size_t t)->bool
  {
    std::pair<double, int> result(1e308, 0);
    for (int k = 0; k < nparam; ++k)
    {
      double *paramPtr = abc_ + k * 3;
      double abc_lm1[4] = { paramPtr[0], paramPtr[1], paramPtr[2], lm1_[i] };
      
      
      // for (int u = 0; u < 4; ++u)
      //   Rcout << abc_lm1[u] << ", ";
      // Rcout << "    ";
      // 
      // 
      // for (int u = 0; u < 3; ++u)
      //   Rcout << abcLB_[u] << ", ";
      // Rcout << "    ";
      // 
      // 
      // for (int u = 0; u < 3; ++u)
      //   Rcout << abcUB_[u] << ", ";
      // Rcout << "    ";
      // 
      // 
      // for (int u = 0; u < distrPtrs_[i].size; ++u)
      //   Rcout << distrPtrs_[i].val[u] << ", ";
      // Rcout << "    ";
      // 
      // 
      // for (int u = 0; u < distrPtrs_[i].size; ++u)
      //   Rcout << distrPtrs_[i].P[u] << ", ";
      // Rcout << "\n\n";
      
      
      auto rst = solveShapeParamsGivenMethodLBFGS<method, RIBlib> (
        abc_lm1, abcLB_, abcUB_, 
        distrPtrs_[i].val, distrPtrs_[i].P, distrPtrs_[i].size, 
        scaleEps, scaleMaxit, opt_[t], max_iterations, hgrad, centralDiff, 
        epsilon, epsilon_rel, m, past, delta, max_submin, max_linesearch, 
        min_step, max_step, ftol, wolfe
      );
      
      
      // Rcout << rst.first << ", " << rst.second << "\n";
      
      
      if (rst.first < result.first) 
      { 
        result = rst;
        std::copy(abc_lm1, abc_lm1 + 4, rstParam_ + i * 4);
      }
      // Rcout << result.second << ", ";
    }
    // Rcout << "\n\n";
    fval_[i] = result.first;
    niter_[i] = result.second;
    return false;
  };
  
  
  cp.parFor(0, fval.size(), runner, grainSize,
            [](std::size_t t)->bool{ return false; },
            [](std::size_t t)->bool{ return false; });
  
  
  return List::create(Named("param") = rstParam, 
                      Named("fval")  = fval, 
                      Named("niter") = niter);
}




// This is suitable for fitting a sequence of PMFs that are similar to
// each other.
template <int method, int RIBlib> 
List LBFGSBtrbFitListGivenMethodSequential(
    NumericMatrix abc,
    NumericVector lm1,
    List empDistrList,
    int beginIndex,
    NumericVector abcLB,
    NumericVector abcUB,
    double scaleEps,
    int scaleMaxit,
    int max_iterations,
    int maxCore,
    double hgrad,
    bool centralDiff,
    int m,
    double epsilon,
    double epsilon_rel,
    int past,
    double delta,
    int max_submin,
    int max_linesearch,
    double min_step,
    double max_step,
    double ftol,
    double wolfe
)
{
  
  if (abc.nrow() != 3) stop("abc's number of rows != 3");
  
  
  if (lm1.size() != empDistrList.size())
    stop("empDistrList's size unequal to lm1's size.");
  double *lm1_ = &lm1[0];
  
  
  const int Ndistr = empDistrList.size();
  const int NparamSet = abc.ncol();
  double *paramPtr = &abc[0];
  
  
  // double *abc_ = &abc[0], *abcLB_ = nullptr, *abcUB_ = nullptr;
  double *abcLB_ = nullptr, *abcUB_ = nullptr;
  if (abcLB.size() != 0) abcLB_ = &abcLB[0];
  if (abcUB.size() != 0) abcUB_ = &abcUB[0];
  
  
  struct DistrPtr { int size; double *val, *P; };
  std::vector<DistrPtr> distrPtrs(Ndistr);
  for (int i = 0; i < Ndistr; ++i)
  { 
    List distr = empDistrList[i];
    NumericVector val = distr[0], P = distr[1];
    distrPtrs[i].size = val.size();
    distrPtrs[i].val = &val[0];
    distrPtrs[i].P = &P[0];
  } 
  DistrPtr *distrPtrs_ = &distrPtrs[0];
  
  
  maxCore = std::min(maxCore, NparamSet);
  CharlieThreadPool cp(std::move(maxCore));
  // maxCore = cp.maxCore;
  std::vector<LBFGSBcharlie> opt(maxCore);
  LBFGSBcharlie *opt_ = &opt[0];
  int grainSize = NparamSet / (maxCore * maxCore * maxCore) + 1;
  
  
  std::vector<double> allFittedParam(NparamSet * int64_t(Ndistr) * 4);
  double *fittedParam = allFittedParam.data();
  double *fittedFval = fittedParam + NparamSet * int64_t(Ndistr) * 3;
  // The 3 parameters fitted for the k-th distribution using 
  //   the i-th initial set starts at 
  //   fittedParam + i * Ndistr * 3 + k * 3 = fittedParam + (i * Ndistr + k) * 3
  // The final function value for fitting the k-th distribution using the 
  //   i-th initial set starts at fittedFval + i * Ndistr + k;
  
  
  std::vector<int> allFitNiter(NparamSet * int64_t(Ndistr));
  int *fitNiter = allFitNiter.data();
  // The iteration for fitting the k-th distribution using the i-th
  // initial set starts at fitNiter + i * Ndistr + k.
  
  
  // i is the i-th set of initial parameters.
  auto runner = [&](std::size_t i, std::size_t t)->bool
  { 
    // double abc_lm1[4] = { 
    //   paramPtr[i * 3], paramPtr[i * 3 + 1], paramPtr[i * 3 + 2], lm1_[0] };
    double abc_lm1[4] = { 
      paramPtr[i * 3], paramPtr[i * 3 + 1], paramPtr[i * 3 + 2], 0.0 };
    
    
    auto fk = [&](int k)->void
    {
      abc_lm1[3] = lm1_[k];
      auto rst = solveShapeParamsGivenMethodLBFGS<method, RIBlib> (
        abc_lm1, abcLB_, abcUB_, 
        distrPtrs_[k].val, distrPtrs_[k].P, distrPtrs_[k].size, 
        scaleEps, scaleMaxit, opt_[t], max_iterations, hgrad, centralDiff, 
        epsilon, epsilon_rel, m, past, delta, max_submin, max_linesearch, 
        min_step, max_step, ftol, wolfe
      );
      auto u = i * Ndistr + k;
      std::copy(abc_lm1, abc_lm1 + 3, fittedParam + u * 3);
      fittedFval[u] = rst.first;
      fitNiter[u] = rst.second;
    };
    
    
    fk(beginIndex);
    double abc_lm1Resv[4];
    std::copy(abc_lm1, abc_lm1 + 3, abc_lm1Resv);
    for (int k = beginIndex + 1; k < Ndistr; ++k) fk(k);
    std::copy(abc_lm1Resv, abc_lm1Resv + 3, abc_lm1); // Must copy back! Inside fk abc_lm1 is used. 
    for (int k = beginIndex - 1; k >= 0; --k) fk(k);
    
    
    /*
    for (int k = beginIndex; k < Ndistr; ++k) fk(k);
    abc_lm1[0] = paramPtr[i * 3]; 
    abc_lm1[1] = paramPtr[i * 3 + 1]; 
    abc_lm1[2] = paramPtr[i * 3 + 2]; 
    for (int k = beginIndex; k >= 0; --k) fk(k);
    */
    
    
    /*
    for (int k = 0; k < Ndistr; ++k)
    {
      abc_lm1[3] = lm1_[k];
      auto rst = solveShapeParamsGivenMethodLBFGS<method, RIBlib> (
        abc_lm1, abcLB_, abcUB_, 
        distrPtrs_[k].val, distrPtrs_[k].P, distrPtrs_[k].size, 
        scaleEps, scaleMaxit, opt_[t], max_iterations, hgrad, centralDiff, 
        epsilon, epsilon_rel, m, past, delta, max_submin, max_linesearch, 
        min_step, max_step, ftol, wolfe
      );
      auto u = i * Ndistr + k;
      std::copy(abc_lm1, abc_lm1 + 3, fittedParam + u * 3);
      fittedFval[u] = rst.first;
      fitNiter[u] = rst.second;
    }
    */
    
    
    return false;
  };  
  
  
  cp.parFor(0, NparamSet, runner, grainSize,
            [](std::size_t t)->bool{ return false; },
            [](std::size_t t)->bool{ return false; });
  
  
  NumericVector fval(Ndistr);
  std::fill(fval.begin(), fval.end(), 1e300);
  IntegerVector niter(Ndistr);
  NumericMatrix rstParam(4, Ndistr);
  double *rstParam_ = &rstParam[0];
  
  
  // The final function value for fitting the k-th distribution using the 
  //   i-th initial set starts at fittedFval + i * Ndistr + k;
  // fittedFval + i * Ndistr + k;
  // int *whichBest = &niter[0];
  for (int i = 0; i < NparamSet; ++i)
  {
    for (int k = 0; k < Ndistr; ++k)
    {
      int64_t u = i * int64_t(Ndistr) + k;
      constexpr const int64_t four = 4;
      if (fittedFval[u] < fval[k])
      {
        fval[k] = fittedFval[u];
        niter[k] = fitNiter[u];
        std::copy(fittedParam + u * 3, fittedParam + u * 3 + 3,
                  rstParam_ + four * k);
      }
    }
  }
  
  
  for (int k = 0; k < Ndistr; ++k)
  {
    const int64_t four = 4;
    rstParam_[four * k + 3] = lm1[k];
  }
  
  
  return List::create(Named("param") = rstParam, 
                      Named("fval")  = fval, 
                      Named("niter") = niter);
} 




template <int method, int RIBlib> 
List LBFGSBtrbFitListGivenMethod(
    NumericMatrix abc,
    NumericVector lm1,
    List empDistrList,
    NumericVector abcLB,
    NumericVector abcUB,
    double scaleEps,
    int scaleMaxit,
    int max_iterations,
    int maxCore,
    int sequentialUpdate, // -1 means no sequential update, otherwise it is the zero-based starting index.
    double hgrad,
    bool centralDiff,
    int m,
    double epsilon,
    double epsilon_rel,
    int past,
    double delta,
    int max_submin,
    int max_linesearch,
    double min_step,
    double max_step,
    double ftol,
    double wolfe
)
{
  if (sequentialUpdate < 0)
    return LBFGSBtrbFitListGivenMethodPlain<method, RIBlib>(
      abc, lm1, empDistrList, abcLB, abcUB, scaleEps, scaleMaxit, max_iterations,
      maxCore, hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta,
      max_submin, max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  return LBFGSBtrbFitListGivenMethodSequential<method, RIBlib>(
    abc, lm1, empDistrList, sequentialUpdate, abcLB, abcUB, scaleEps, 
    scaleMaxit, max_iterations, maxCore, hgrad, centralDiff, m, epsilon, 
    epsilon_rel, past, delta,
    max_submin, max_linesearch, min_step, max_step, ftol, wolfe);
}








template <int RIBlib>
List LBFGSBtrbFitList(
    NumericMatrix abc,
    NumericVector lm1,
    List empDistrList,
    NumericVector abcLB,
    NumericVector abcUB,
    double scaleEps,
    int scaleMaxit,
    String distanceFun,
    int maxCore,
    int sequentialUpdate,
    double hgrad,
    bool centralDiff,
    int max_iterations,
    int m,
    double epsilon,
    double epsilon_rel,
    int past,
    double delta,
    int max_submin,
    int max_linesearch,
    double min_step,
    double max_step,
    double ftol,
    double wolfe
)
{
  
  List rst = List::create();
  
  
  if (distanceFun == "likelihood") 
    rst = LBFGSBtrbFitListGivenMethod<0, RIBlib> (
    abc, lm1, empDistrList, abcLB, abcUB, scaleEps, scaleMaxit, max_iterations,
    maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta,
    max_submin, max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  else if (distanceFun == "likelihoodDiscrete") 
    rst = LBFGSBtrbFitListGivenMethod<1, RIBlib> (
    abc, lm1, empDistrList, abcLB, abcUB, scaleEps, scaleMaxit, max_iterations,
    maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta,
    max_submin, max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  else if (distanceFun == "Komogorov") 
    rst = LBFGSBtrbFitListGivenMethod<2, RIBlib> (
    abc, lm1, empDistrList, abcLB, abcUB, scaleEps, scaleMaxit, max_iterations,
    maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta,
    max_submin, max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  else if (distanceFun == "EuclideanCDF") 
    rst = LBFGSBtrbFitListGivenMethod<3, RIBlib> (
    abc, lm1, empDistrList, abcLB, abcUB, scaleEps, scaleMaxit, max_iterations,
    maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta,
    max_submin, max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  else if (distanceFun == "unbiasedLikelihoodDiscrete") 
    rst = LBFGSBtrbFitListGivenMethod<4, RIBlib> (
      abc, lm1, empDistrList, abcLB, abcUB, scaleEps, scaleMaxit, max_iterations,
      maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta,
      max_submin, max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  
  else if (distanceFun == "unbiasedKomogorov") 
    rst = LBFGSBtrbFitListGivenMethod<5, RIBlib> (
      abc, lm1, empDistrList, abcLB, abcUB, scaleEps, scaleMaxit, max_iterations,
      maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta,
      max_submin, max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  
  else if (distanceFun == "unbiasedEuclideanCDF") 
    rst = LBFGSBtrbFitListGivenMethod<6, RIBlib> (
      abc, lm1, empDistrList, abcLB, abcUB, scaleEps, scaleMaxit, max_iterations,
      maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta,
      max_submin, max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  
  else stop("Distance function not implemented");
  
  
  return rst;
}




//' Fit TrBs
//'
//' Fit TrBs to a list of PMFs using \href{https://en.wikipedia.org/wiki/Limited-memory_BFGS#L-BFGS-B}{L-BFGS-B}. 
//' 
//' @param abc  Initializations of \code{a,b,c} as a \eqn{3\times K} numeric 
//' matrix. Here \eqn{K} is the number of initializations. The function will 
//' fit each distribution using the \eqn{K} initializations, and then return  
//' the optimized parameters that lead to the best 
//' objective function value.
//' 
//' @param lm1  A numeric vector of target limited means for constraining 
//' the TrBs.
//' 
//' @param empDistrList  A list of empirical PMFs. Each PMF is a list of two 
//' numeric vectors, the support and the probabilities.
//' 
//' @param abcLB  A numeric vector as the lower bounds on \code{a,b,c}.
//' An emptry vector \code{numeric(0)} implies no lower bounds. 
//' Supplying the lower bounds is strongly recommended. 
//' Default to \code{c(1.01,0.1,0.1)}. Note that \code{a<=1} implies 
//' a TrB distribution with infinite mean.
//' 
//' @param abcUB  A numeric vector as the upper bounds on \code{a,b,c}. 
//' An emptry vector \code{numeric(0)} imples no upper bounds. 
//' Supplying the upper bounds is strongly recommended. 
//' Default to \code{c(30,30,30)}.
//' 
//' @param scaleEps  Stop Newton's solver for parameter \code{d} if 
//' the difference between TrB's current limited mean and
//' its target limited mean becomes less than \code{scaleEps}.
//' Default to \code{1e-8}.
//' 
//' @param scaleMaxit  Stop Newton's solver for parameter \code{d} if 
//' the number of iterations has reached \code{scaleMaxit}. 
//' Default to \code{100}.
//' 
//' @param distanceFun  A string. Name of the distance function for optimization.
//' \itemize{
//' 
//' \item{\code{"likelihood"}:}{ Negative log-likelihood.}
//' 
//' \item{\code{"likelihoodDiscrete"}:}{ Discretize the TrB and then compute 
//' the negative log-likelihood, aka cross-entropy. Much slower due to 
//' frequent calls to TrB's CDF. No improvement over 
//' \code{"likelihood"} was found.}
//'   
//' \item{\code{"Komogorov"}:}{\href{https://en.wikipedia.org/wiki/Smooth_maximum}{ Smooth maximum} 
//' distance between CDFs. Slow. The smooth maximum ensures differentiability. }
//' 
//' \item{\code{"EuclideanCDF"}:}{ Euclidean distance between CDFs. Slow.}
//' 
//' \item{\code{"unbiasedLikelihoodDiscrete"}:}{ A research effort misled by the 
//'   "unbiased" discretization method described in R package \code{actuar}. 
//'   Such discretization preserves the PDF's mean, but the probabilities 
//'   will not sum up to 1. Use with caution.}
//'   
//' \item{\code{"unbiasedKomogorov"}:}{ A research effort misled by the 
//'   "unbiased" discretization method described in R package \code{actuar}. 
//'   Such discretization preserves the PDF's mean, but the probabilities 
//'   will not sum up to 1. Use with caution.}
//'   
//' \item{\code{"unbiasedEuclideanCDF"}:}{ A research effort misled by the 
//'   "unbiased" discretization method described in R package \code{actuar}. 
//'   Such discretization preserves the PDF's mean, but the probabilities 
//'   will not sum up to 1. Use with caution.}
//' }
//' 
//' 
//' @param max_iterations  An L-BFGS-B algorithm parameter. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @param maxCore  The number of CPU cores for use. Default to 1000. A number
//'   greater than the actual number of logical processors will be capped.
//'   The bidirectional sequential fitting, which is the most likely to be used, 
//'   is single-threaded if the initialization \code{abc} is singular --- when
//'   \eqn{K=1}.
//'   
//' @param RIBlib  A string to specify which library should be used to 
//' compute the \href{https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function}{Regularized Incomplete Beta function}
//' that TrB CDF depends on. 
//' \itemize{
//' \item{\code{"Numerical Recipes"}}{ implements the textbook algorithm in
//' \href{https://e-maxx.ru/bookz/files/numerical_recipes.pdf}{Chapter 6.4, Numerical Recipes}.
//' \href{https://live.boost.org/doc/libs/1_54_0/libs/math/doc/html/math_toolkit/sf_beta/ibeta_function.html}{Boost}, 
//' which claims the textbook algorithm's series expansion converges less 
//' rapidly than its own approach, actually runs considerably slower because 
//' it pursues much higher numeric precision.
//' Such extreme precision turns out largely unncessary in our setting. 
//' \code{"Numerical Recipes"} is the default and 
//' the recommended choice.}
//' 
//' \item{\code{"R::pbeta"}}{ calls R's Beta CDF. It is slower than \code{"Numerical Recipes"}
//' but faster than \code{"Boost"}. In numerically pathological regions, 
//' \code{"R::pbeta"} is more precise than \code{"Numerical Recipes"} 
//' but less so than \code{"Boost"}.}
//' 
//' \item{\code{"Boost"}}{ calls \href{https://live.boost.org/doc/libs/1_54_0/libs/math/doc/html/math_toolkit/sf_beta/ibeta_function.html}{boost::math::ibeta}.
//' It is the slowest but yields much higher precision especially 
//' in numerically pathological regions.}
//' }
//' 
//' @param sequentialUpdate  One-based index of the first PMF to be fitted 
//' during bi-directional sequential fitting. \code{-1} avoids using 
//' the bi-directional sequential fitting.
//' 
//' @param hgrad  A nonnegative number as the \eqn{\Delta} for finite differencing. 
//' Default to \code{0} which implies to
//' use the internal default values: \code{6.055e-6} for central difference, 
//' \code{1.49e-8} for forward difference. The numbers are 
//' cubic and square roots of machine precision respectively.
//' 
//' @param centralDiff  \code{TRUE} to use central difference for finite 
//' differencing. 
//' \code{FALSE} to use forward difference. Default to \code{TRUE} 
//' which is also strongly recommended.
//' 
//' @param m  Parameter in L-BFGS-B. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @param epsilon  Parameter in L-BFGS-B. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @param epsilon_rel  Parameter in L-BFGS-B. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @param past  Parameter in L-BFGS-B. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @param delta Parameter in L-BFGS-B. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @param max_submin Parameter in L-BFGS-B. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @param max_linesearch Parameter in L-BFGS-B. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @param min_step Parameter in L-BFGS-B. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @param max_step  Parameter in L-BFGS-B. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @param ftol Parameter in L-BFGS-B. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @param wolfe Parameter in L-BFGS-B. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @return A list of 3 objects:
//' \itemize{
//' \item{\code{param}:}{ A numeric matrix of 4 rows. \code{param[, i]} 
//' contains the optimized \code{a,b,c} and the input \code{lm1} for 
//' the \code{i}-th empirical PMF.}
//' 
//' \item{\code{fval}:}{ The final objective function values. A value of 
//' \code{1e300} implies the optimization has not converged and the 
//' corresponding parameters should be reanalyzed or discarded.}
//' 
//' \item{\code{niter}:}{ The number of iterations until convergence.}
//' }
//' 
//' 
//' @inherit  computeWindows details
//' 
//' @example  inst/examples/LBFGSBtrbFitList.R
//' 
//' 
// [[Rcpp::export]]
List LBFGSBtrbFitList(
    NumericMatrix abc,
    NumericVector lm1,
    List empDistrList,
    NumericVector abcLB = NumericVector::create(1.01, 0.1, 0.1),
    NumericVector abcUB = NumericVector::create(30, 30, 30),
    double scaleEps = 1e-8, 
    int scaleMaxit = 100,
    String distanceFun = "likelihood",
    int max_iterations = 100,
    int maxCore = 15,
    String RIBlib = "Numerical Recipes",
    int sequentialUpdate = -1,
    double hgrad = 0,
    bool centralDiff = true,
    int m = 6,
    double epsilon = 1e-5,
    double epsilon_rel = 1e-5,
    int past = 1,
    double delta = 1e-10,
    int max_submin = 10,
    int max_linesearch = 20,
    double min_step = 1e-20,
    double max_step = 1e+20,
    double ftol = 1e-4,
    double wolfe = 0.9
)
{
  
  // If sequential update, sequentialUpdate would be the 1-based starting index 
  //   for bi-directional sequential update. 
  if (sequentialUpdate > 0) sequentialUpdate -= 1;
  
  
  if (RIBlib == "Boost") return LBFGSBtrbFitList<0> (
    abc,lm1,empDistrList, abcLB, abcUB, scaleEps,scaleMaxit,distanceFun,
    max_iterations, maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, 
    epsilon_rel, past, delta, max_submin, max_linesearch, min_step, 
    max_step, ftol, wolfe);
  
  
  if (RIBlib == "R::pbeta") return LBFGSBtrbFitList<1> (
    abc,lm1,empDistrList, abcLB,abcUB,scaleEps,scaleMaxit,distanceFun,
    max_iterations, maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, 
    epsilon_rel, past, delta, max_submin, max_linesearch, min_step, 
    max_step, ftol, wolfe);
  
  
  if (RIBlib == "Numerical Recipes") return LBFGSBtrbFitList<2> (
    abc,lm1,empDistrList, abcLB,abcUB,scaleEps,scaleMaxit,distanceFun,
    max_iterations, maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, 
    epsilon_rel, past, delta, max_submin, max_linesearch, min_step, 
    max_step, ftol, wolfe);
  
  
  stop("Method of computing Regularized Incomplete Beta function is not implemented.");
  
  
  return List::create();
}









#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// LBFGSBtrbFit
List LBFGSBtrbFit(NumericVector abc_lm1, DataFrame empDistr, NumericVector abcLB, NumericVector abcUB, double scaleEps, int scaleMaxit, String distanceFun, int max_iterations, String RIBlib, double hgrad, bool centralDiff, int m, double epsilon, double epsilon_rel, int past, double delta, int max_submin, int max_linesearch, double min_step, double max_step, double ftol, double wolfe);
RcppExport SEXP sourceCpp_32_LBFGSBtrbFit(SEXP abc_lm1SEXP, SEXP empDistrSEXP, SEXP abcLBSEXP, SEXP abcUBSEXP, SEXP scaleEpsSEXP, SEXP scaleMaxitSEXP, SEXP distanceFunSEXP, SEXP max_iterationsSEXP, SEXP RIBlibSEXP, SEXP hgradSEXP, SEXP centralDiffSEXP, SEXP mSEXP, SEXP epsilonSEXP, SEXP epsilon_relSEXP, SEXP pastSEXP, SEXP deltaSEXP, SEXP max_subminSEXP, SEXP max_linesearchSEXP, SEXP min_stepSEXP, SEXP max_stepSEXP, SEXP ftolSEXP, SEXP wolfeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type abc_lm1(abc_lm1SEXP);
    Rcpp::traits::input_parameter< DataFrame >::type empDistr(empDistrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type abcLB(abcLBSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type abcUB(abcUBSEXP);
    Rcpp::traits::input_parameter< double >::type scaleEps(scaleEpsSEXP);
    Rcpp::traits::input_parameter< int >::type scaleMaxit(scaleMaxitSEXP);
    Rcpp::traits::input_parameter< String >::type distanceFun(distanceFunSEXP);
    Rcpp::traits::input_parameter< int >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< String >::type RIBlib(RIBlibSEXP);
    Rcpp::traits::input_parameter< double >::type hgrad(hgradSEXP);
    Rcpp::traits::input_parameter< bool >::type centralDiff(centralDiffSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon_rel(epsilon_relSEXP);
    Rcpp::traits::input_parameter< int >::type past(pastSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type max_submin(max_subminSEXP);
    Rcpp::traits::input_parameter< int >::type max_linesearch(max_linesearchSEXP);
    Rcpp::traits::input_parameter< double >::type min_step(min_stepSEXP);
    Rcpp::traits::input_parameter< double >::type max_step(max_stepSEXP);
    Rcpp::traits::input_parameter< double >::type ftol(ftolSEXP);
    Rcpp::traits::input_parameter< double >::type wolfe(wolfeSEXP);
    rcpp_result_gen = Rcpp::wrap(LBFGSBtrbFit(abc_lm1, empDistr, abcLB, abcUB, scaleEps, scaleMaxit, distanceFun, max_iterations, RIBlib, hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta, max_submin, max_linesearch, min_step, max_step, ftol, wolfe));
    return rcpp_result_gen;
END_RCPP
}
// LBFGSBtrbFitList
List LBFGSBtrbFitList(NumericMatrix abc, NumericVector lm1, List empDistrList, NumericVector abcLB, NumericVector abcUB, double scaleEps, int scaleMaxit, String distanceFun, int max_iterations, int maxCore, String RIBlib, int sequentialUpdate, double hgrad, bool centralDiff, int m, double epsilon, double epsilon_rel, int past, double delta, int max_submin, int max_linesearch, double min_step, double max_step, double ftol, double wolfe);
RcppExport SEXP sourceCpp_32_LBFGSBtrbFitList(SEXP abcSEXP, SEXP lm1SEXP, SEXP empDistrListSEXP, SEXP abcLBSEXP, SEXP abcUBSEXP, SEXP scaleEpsSEXP, SEXP scaleMaxitSEXP, SEXP distanceFunSEXP, SEXP max_iterationsSEXP, SEXP maxCoreSEXP, SEXP RIBlibSEXP, SEXP sequentialUpdateSEXP, SEXP hgradSEXP, SEXP centralDiffSEXP, SEXP mSEXP, SEXP epsilonSEXP, SEXP epsilon_relSEXP, SEXP pastSEXP, SEXP deltaSEXP, SEXP max_subminSEXP, SEXP max_linesearchSEXP, SEXP min_stepSEXP, SEXP max_stepSEXP, SEXP ftolSEXP, SEXP wolfeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type abc(abcSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lm1(lm1SEXP);
    Rcpp::traits::input_parameter< List >::type empDistrList(empDistrListSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type abcLB(abcLBSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type abcUB(abcUBSEXP);
    Rcpp::traits::input_parameter< double >::type scaleEps(scaleEpsSEXP);
    Rcpp::traits::input_parameter< int >::type scaleMaxit(scaleMaxitSEXP);
    Rcpp::traits::input_parameter< String >::type distanceFun(distanceFunSEXP);
    Rcpp::traits::input_parameter< int >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type maxCore(maxCoreSEXP);
    Rcpp::traits::input_parameter< String >::type RIBlib(RIBlibSEXP);
    Rcpp::traits::input_parameter< int >::type sequentialUpdate(sequentialUpdateSEXP);
    Rcpp::traits::input_parameter< double >::type hgrad(hgradSEXP);
    Rcpp::traits::input_parameter< bool >::type centralDiff(centralDiffSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon_rel(epsilon_relSEXP);
    Rcpp::traits::input_parameter< int >::type past(pastSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type max_submin(max_subminSEXP);
    Rcpp::traits::input_parameter< int >::type max_linesearch(max_linesearchSEXP);
    Rcpp::traits::input_parameter< double >::type min_step(min_stepSEXP);
    Rcpp::traits::input_parameter< double >::type max_step(max_stepSEXP);
    Rcpp::traits::input_parameter< double >::type ftol(ftolSEXP);
    Rcpp::traits::input_parameter< double >::type wolfe(wolfeSEXP);
    rcpp_result_gen = Rcpp::wrap(LBFGSBtrbFitList(abc, lm1, empDistrList, abcLB, abcUB, scaleEps, scaleMaxit, distanceFun, max_iterations, maxCore, RIBlib, sequentialUpdate, hgrad, centralDiff, m, epsilon, epsilon_rel, past, delta, max_submin, max_linesearch, min_step, max_step, ftol, wolfe));
    return rcpp_result_gen;
END_RCPP
}
