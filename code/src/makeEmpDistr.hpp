
// [[Rcpp::export]]
DataFrame makeEmpDistr(
    NumericVector x, DataFrame pmf, double w, 
    int rstSize, String regridMethod = "r4", 
    double fixedMin = 1e300, double fixedMax = 1e300,
    double biasCorrectionMultiplier = 1.0)
{
  
  MergeRegrid mr;
  NumericVector y = pmf[0], yp = pmf[1];
  NumericVector rst(rstSize), rstP(rstSize);
  if (regridMethod == "lr")
  {
    mr.operator()<0, int, double, double, double, double, double, double>(
        &x[0], x.size(), &y[0], &yp[0], y.size(),
        &rst[0], &rstP[0], rstSize, w, fixedMin, fixedMax, biasCorrectionMultiplier);
  }
  else if (regridMethod == "lmm")
  {
    mr.operator()<1, int, double, double, double, double, double, double>(
        &x[0], x.size(), &y[0], &yp[0], y.size(),
        &rst[0], &rstP[0], rstSize, w, fixedMin, fixedMax, biasCorrectionMultiplier);
  }
  else if (regridMethod == "r4")
  {
    mr.operator()<2, int, double, double, double, double, double, double>(
        &x[0], x.size(), &y[0], &yp[0], y.size(),
        &rst[0], &rstP[0], rstSize, w, fixedMin, fixedMax, biasCorrectionMultiplier);
  }
  else stop("Regrid method not implemented.");
  
  
  return List::create(Named("val") = rst, Named("P") = rstP);
}




// window is a 2-row integer matrix. The first row is the 1-based index
// of the first element in the window over X. The second row is the 1-based 
// index of the last element in the window over X.
// w is the weight on PMFs.
template <int regridMethod>
List makeEmpDistrList(NumericVector X, 
                      IntegerMatrix windows,
                      IntegerVector rstSizes, 
                      int maxCore, double fixedMin, double fixedMax,
                      NumericVector biasCorrectionMultiplier)
{
  const int Nwindow = windows.ncol();
  if (rstSizes.size() != 1 and rstSizes.size() != Nwindow) 
    stop("rstSizes.size() != windows.ncol().");
  
  
  for (int u: rstSizes) 
  { 
    if (u < 2) stop("rstSizes has element less than 2."); 
  }
  
  
  int *win = &windows[0];
  double *x = &X[0];
  CharlieThreadPool cp(maxCore);
  maxCore = cp.maxCore;
  std::vector<MergeRegrid> MR(maxCore);
  List rst(Nwindow);
  int rstSze = rstSizes.size();
  
  
  std::vector<double*> ptrs(Nwindow * 2);
  double **valPtrs = ptrs.data();
  double **Pptrs = valPtrs + Nwindow;
  int *distrSizes = &rstSizes[0];
  for (int i = 0; i < Nwindow; ++i)
  {
    NumericVector val(distrSizes[i % rstSze]);
    valPtrs[i] = &val[0];
    NumericVector P(distrSizes[i % rstSze]);
    Pptrs[i] = &P[0];
    rst[i] = List::create(Named("val") = val, Named("P") = P);
  }
  
  
  auto f = [&](
    std::size_t i, std::size_t t)->bool
  {
    MergeRegrid &mr = MR[t];
    double *xbegin = win[2 * i] - 1 + x;
    int xsize = win[2 * i + 1] - win[2 * i] + 1;
    int bmsize = biasCorrectionMultiplier.size();
    mr.operator() <regridMethod, int, double, double, double, double, double, double> (
        xbegin, xsize, nullptr, nullptr, 0, 
        valPtrs[i], Pptrs[i], distrSizes[i % rstSze], 0.0, fixedMin, fixedMax,
        biasCorrectionMultiplier[i % bmsize]
    );
    return false;
  };
  
  
  int grainSize = Nwindow / (maxCore * maxCore) + 1;
  cp.parFor(0, Nwindow, f, grainSize, 
            [](std::size_t t)->bool { return false; },
            [](std::size_t t)->bool { return false; });
  
  
  return rst;
}




// [[Rcpp::export]]
List makeEmpDistrList(NumericVector X, 
                      IntegerMatrix windows,
                      IntegerVector rstSizes, 
                      String regridMethod = "r4",
                      int maxCore = 1000, 
                      double fixedMin = 1e300,
                      double fixedMax = 1e300,
                      NumericVector biasCorrectionMultiplier = 1.0)
{
  
  if (regridMethod == "lr")
    return makeEmpDistrList<0>(X, windows, rstSizes, maxCore, fixedMin, fixedMax,
                               biasCorrectionMultiplier);
  
  
  if (regridMethod == "lmm")
    return makeEmpDistrList<1>(X, windows, rstSizes, maxCore, fixedMin, fixedMax,
                               biasCorrectionMultiplier);
  
  
  if (regridMethod == "r4")
    return makeEmpDistrList<2>(X, windows, rstSizes, maxCore, fixedMin, fixedMax,
                               biasCorrectionMultiplier);
  
  
  stop("Regrid method not implemented");
  return List::create();
}
  









