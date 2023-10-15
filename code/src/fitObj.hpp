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
  CharlieThreadPool cp(maxCore);
  maxCore = cp.maxCore;
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
  CharlieThreadPool cp(maxCore);
  maxCore = cp.maxCore;
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
    int max_iterations,
    int maxCore,
    int sequentialUpdate,
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
//' @param abc  Initializations of \code{a, b, c} as a numeric matrix. It can
//' be a \eqn{3 \times 1} matrix, or a \eqn{3 \times }\code{length(lm1)} matrix.
//' Having only 1 column implies equal initializations for all TrB fitting.
//' 
//' @param lm1  A numeric vector of target limited means for constraining the TrBs.
//' 
//' @param empDistrList  A list of empirical PMFs. Each PMF is a list of two 
//' numeric vectors, the support and the probabilities.
//' 
//' @param abcLB  A numeric vector of the lower bounds on \code{a, b, c}.
//' Default to \code{numeric(0)} which implies \eqn{-\infty}. Supplying the lower 
//' bounds is strongly recommended.
//' 
//' @param abcUB  A numeric vector of the upper bounds on \code{a, b, c}. 
//' Default to \code{numeric(0)} which imples \eqn{\infty}. Supplying the upper
//' bounds is strongly recommended.
//' 
//' @param scaleEps  If the difference between TrB's actual limited mean and
//' its target limited mean is less than \code{scaleEps}, stop the Newton's
//' iteration for solving parameter \code{d}. Default to \code{1e-8}.
//' 
//' @param scaleMaxit  If the number of Newton's iterations for solving \code{d}
//' reaches \code{scaleMaxit}, stop. Default to \code{100}.
//' 
//' @param distanceFun  Name of the distance function for optimization.
//' \itemize {
//' 
//' \item\code{"likelihood"}: Negative log-likelihood.
//' 
//' \item\code{"likelihoodDiscrete"}: Discretize the TrB and then compute the negative 
//'   log-likelihood, aka cross-entropy. Much slower due to frequent calls to
//'   the TrB's CDF. No improvements over "likelihood" were found.
//'   
//' \item\code{"Komogorov"}: Smooth maximum distance between CDFs. Slow.
//' 
//' \item\code{"EuclideanCDF"}: Euclidean distance between CDFs. Slow.
//' 
//' \item\code{"unbiasedLikelihoodDiscrete"}: A research effort misled by the 
//'   "unbiased" discretization method described in R package \code{actuar}. 
//'   Such discretization preserves the PDF's mean, but the probabilities 
//'   will not sum up to 1. Not for use.
//'   
//' \item\code{"unbiasedKomogorov"}: A research effort misled by the 
//'   "unbiased" discretization method described in R package \code{actuar}. 
//'   Such discretization preserves the PDF's mean, but the probabilities 
//'   will not sum up to 1. Not for use.
//'   
//' \item\code{"unbiasedEuclideanCDF"}: A research effort misled by the 
//'   "unbiased" discretization method described in R package \code{actuar}. 
//'   Such discretization preserves the PDF's mean, but the probabilities 
//'   will not sum up to 1. Not for use.
//' }
//' 
//' 
//' @param max_iterations  A L-BFGS-B algorithm parameter. Do not change the
//' default value before consulting \href{https://github.com/yixuan/LBFGSpp}{LBFGSpp}.
//' 
//' @param maxCore  The number of CPU cores for use. Default to 1000. A number
//'   greater than the actual number of logical processors will be capped.
//'   The bidirectional sequential fitting, which is the most likely to be used, 
//'   is single-threaded however.
//'   
//' @param RIBlib  Use what library to compute the Regularized Incomplete Beta
//' function which TrB CDF depends on. 
//' \itemize{
//' \item\code{"Numerical Recipes"} implements the textbook algorithm in
//' "Numerical Recipes". \href{https://live.boost.org/doc/libs/1_54_0/libs/math/doc/html/math_toolkit/sf_beta/ibeta_function.html}{Boost},
//' claims the textbook's series expansion converges slower, but Boost's 
//' runtime takes much longer as it pursues extremely high precision, which is 
//' unnecessary for our problem. \code{"Numerical Recipes"} is the default and 
//' the recommended choice.
//' 
//' \item\code{"R::pbeta"} calls R's Beta CDF. It is slower than \code{"Numerical Recipes"}
//' but faster than \code{"Boost"}. In numerically pathological region, \code{"R::pbeta"}'s
//' precision is higher than \code{"Numerical Recipes"} but lower than \code{"Boost"}.
//' 
//' \item\code{"Boost"} calls \href{https://live.boost.org/doc/libs/1_54_0/libs/math/doc/html/math_toolkit/sf_beta/ibeta_function.html}{Incomplete Beta Functions}.
//' It is the slowest but has the highest precision among the three libraries.
//' }
//' 
//' @param sequentialUpdate  One-based index of the first PMF to be fitted during bi-directional
//' sequential fitting. \code{-1} implies sequential fitting should not be used.
//' 
//' @param hgrad  \eqn{\Delta} for finite difference. \code{0} implies using the
//' internal default values.
//' 
//' @param centralDiff  TRUE if central difference should be used. FALSE implies
//' forward difference.
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
//' @return A list of three:
//' \itemize{
//' \item\code{$param}: A numeric matrix of four rows. \code{param[, i]} is the
//' optimized \code{a, b, c} and \code{lm1} for the \code{i-th} empirical PMF.
//' 
//' \item\code{$fval}: The final objective function values. A value of 1e300
//' implies the optimization is not converged and the corresponding parameters
//' should be reconsidered or discarded.
//' 
//' \item\code{$niter}: The number of iterations for L-BFGS-B to converge.
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
    NumericVector abcLB = NumericVector(0),
    NumericVector abcUB = NumericVector(0),
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
    max_iterations, maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, epsilon_rel,
    past, delta, max_submin, max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  if (RIBlib == "R::pbeta") return LBFGSBtrbFitList<1> (
    abc,lm1,empDistrList, abcLB,abcUB,scaleEps,scaleMaxit,distanceFun,
    max_iterations, maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, epsilon_rel,
    past, delta, max_submin, max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  if (RIBlib == "Numerical Recipes") return LBFGSBtrbFitList<2> (
    abc,lm1,empDistrList, abcLB,abcUB,scaleEps,scaleMaxit,distanceFun,
    max_iterations, maxCore, sequentialUpdate, hgrad, centralDiff, m, epsilon, epsilon_rel,
    past, delta, max_submin, max_linesearch, min_step, max_step, ftol, wolfe);
  
  
  stop("Method of computing Regularized Incomplete Beta function is not implemented.");
  
  
  return List::create();
}







