

template <int regridMethod>
List mixDistrList(List X, List Y, IntegerVector rstSizes, 
                  NumericVector Yweights, int maxCore)
{
  if (Y.size() != 0 and X.size() != Y.size()) 
    stop("Y.size() != 0 and X.size() != Y.size().");
  for (int u: rstSizes) 
  { 
    if (u < 2) stop("rstSizes has element less than 2."); 
  }
  
  
  CharlieThreadPool cp(maxCore);
  maxCore = cp.maxCore;
  std::vector<MixRegrid> MR(maxCore);
  List rst(X.size());
  int Ndistr = X.size();
  int rstSz = rstSizes.size();
  int wsz = Yweights.size();
  double *yweights = &Yweights[0];
  
  
  std::vector<double*> ptrs(Ndistr * 6);
  double **valPtrs = ptrs.data();
  double **Pptrs = valPtrs + Ndistr;
  int *distrSizes = &rstSizes[0];
  for (int i = 0; i < Ndistr; ++i)
  {
    NumericVector val(distrSizes[i % rstSz]);
    valPtrs[i] = &val[0];
    NumericVector P(distrSizes[i % rstSz]);
    Pptrs[i] = &P[0];
    rst[i] = List::create(Named("val") = val, Named("P") = P);
  }
  
  
  double **xvals = Pptrs + Ndistr, **xps = xvals + Ndistr,
    **yvals = xps + Ndistr, **yps = yvals + Ndistr;
  std::vector<int> xysizes(Ndistr * 2);
  int *Xsizes = xysizes.data(), *Ysizes = Xsizes + Ndistr;
  for (int i = 0; i < Ndistr; ++i)
  {
    List tmp = X[i];
    NumericVector tmp0 = tmp[0]; xvals[i] = &tmp0[0];
    NumericVector tmp1 = tmp[1]; xps  [i] = &tmp1[0];
    Xsizes[i] = tmp0.size();
    List smp = Y[i];
    NumericVector smp0 = smp[0]; yvals[i] = &smp0[0];
    NumericVector smp1 = smp[1]; yps  [i] = &smp1[0];
    Ysizes[i] = smp0.size();
  }
  
  
  auto f = [&](std::size_t i, std::size_t t)->bool
  {
    MixRegrid &mr = MR[t];
    mr.operator()<regridMethod, int, double, double, double, double, double, double> (
          xvals[i], xps[i], Xsizes[i], yvals[i], yps[i], Ysizes[i], 
          valPtrs[i], Pptrs[i], distrSizes[i % rstSz], yweights[i % wsz]);
    return false;
  };
  
  
  int grainSize = Ndistr / (maxCore * maxCore) + 1;
  cp.parFor(0, Ndistr, f, grainSize, 
            [](std::size_t t)->bool { return false; },
            [](std::size_t t)->bool { return false; });
  
  
  return rst;
}




// [[Rcpp::export]]
List mixDistrList(List X, List Y, IntegerVector rstSizes, 
                  NumericVector Yweights,
                  String regridMethod = "r4",
                  int maxCore = 1000)
{
  
  if (regridMethod == "lr")
    return mixDistrList<0>(X, Y, rstSizes, Yweights, maxCore);
  
  
  if (regridMethod == "lmm")
    return mixDistrList<1>(X, Y, rstSizes, Yweights, maxCore);
  
  
  if (regridMethod == "r4")
    return mixDistrList<2>(X, Y, rstSizes, Yweights, maxCore);
  
  
  stop("Regrid method not implemented");
  return List::create();
}




// [[Rcpp::export]]
List mixDistr(List X, List Y, int rstSize, 
              double yweight,
              String regridMethod = "r4")
{
  List Xlist = List::create(X);
  List Ylist = List::create(Y);
  IntegerVector rstSizes = {rstSize};
  NumericVector Yweights = {yweight};
  return mixDistrList(
    Xlist, Ylist, rstSizes, Yweights, regridMethod, 1);
}













