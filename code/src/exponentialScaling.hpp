

struct ExpScaling
{
  Regrid regrider;
  std::vector<double> cntr;
  
  
  template <int regridMethod>
  void operator()(double *val, double *P, int size, 
                double targetMDR, double eps, 
                double *rstval, double *rstP, int maxit)
  {
    cntr.resize(size);
    double *logval = this->cntr.data();
    for (int i = 0; i < size; ++i) 
    {
      if (val[i] <= 0) logval[i] = -1e300;
      else logval[i] = std::log(val[i]); 
    }
    auto f = [this, &P, &size, &targetMDR, &logval](
      double pw)->std::pair<double, double>
      {
        std::pair<double, double> rst(0, 0);
        double &mean = rst.first, &meanDeriv = rst.second;
        for (int i = 0; i < size; ++i)
        {
          if (logval[i] <= -1e300) continue;
          double tpi = std::exp(logval[i] * pw) * P[i];
          mean += tpi;
          meanDeriv += tpi * logval[i];
        }
        rst.first = targetMDR - rst.first; // All vals are in [0, 1], thus decreasing function.
        // Rcout << "meanDeriv = " << meanDeriv << "\n";
        rst.second = -rst.second;
        return rst;
      };
    double pw = newtonRoot(1.0, f, eps, maxit);
    double *valnew = logval;
    for (int i = 0; i < size; ++i) 
      // valnew[i] = std::exp(logval[i] * pw);
      valnew[i] = std::pow(val[i] , pw);
    // Rcout << "pw = " << pw  << "\n";
    makeRegularGrid(valnew[0], valnew[size - 1], rstval, size);
    regrider.operator()<
      regridMethod, int, double, double, 
      double, double, false, true> (
          valnew, P, size, rstval, rstP, size);
  }
};




// This is used to extrapolate the empirical distribution for MDRs below
//   the minimum MDR in the given data.
template <int regridMethod>
List expScalingForMeanMatching(List distlist, NumericVector mdrs, 
                               int maxCore, double eps, int maxIter)
{
  struct dist {int size; double *val, *P;};
  int Ndist = mdrs.size();
  if (distlist.size() != 1 and Ndist != distlist.size()) 
    stop("distlist.size() != mdrs.size()");
  std::vector<dist> cntr(Ndist * 2);
  dist *D = cntr.data(), *Drst = D + Ndist;
  List rst(Ndist);
  int distlistLen = distlist.size();
  for (int i = 0; i < Ndist; ++i)
  {
    List d = distlist[i % distlistLen];
    NumericVector val = d[0], p = d[1];
    D[i].size = val.size();
    D[i].val = &val[0];
    D[i].P = &p[0];
    NumericVector val_(val.size()), P_(val.size());
    Drst[i].size = val_.size();
    Drst[i].val = &val_[0];
    Drst[i].P = &P_[0];
    rst[i] = List::create(Named("val") = val_, Named("P") = P_);
  }
  
  
  CharlieThreadPool cp(maxCore);
  maxCore = cp.maxCore;
  std::vector<ExpScaling> rxp(maxCore);
  auto f = [&](std::size_t i, std::size_t t)->bool
  {
    rxp[t].operator()<regridMethod> (D[i].val, D[i].P, D[i].size,
                    mdrs[i], eps, Drst[i].val, Drst[i].P, maxIter);
    return false;
  };
  int grainSize = Ndist / (maxCore * maxCore * maxCore) + 1;
  cp.parFor(0, Ndist, f, grainSize, 
            [](std::size_t t)->bool { return false; },
            [](std::size_t t)->bool { return false; });
  return rst;
}




// [[Rcpp::export]]
List expScalingForMeanMatching(List distlist, NumericVector mdrs, 
                               String regridMethod = "lr",
                               int maxCore = 1000, double eps = 1e-10, 
                               int maxIter = 100)
{
  if (regridMethod == "lr")
    return expScalingForMeanMatching<0>(
      distlist, mdrs, maxCore, eps, maxIter);
  if (regridMethod == "lmm")
    return expScalingForMeanMatching<1>(
      distlist, mdrs, maxCore, eps, maxIter);
  if (regridMethod == "r4")
    return expScalingForMeanMatching<2>(
      distlist, mdrs, maxCore, eps, maxIter);
  stop("Regrid method not implemented.");
  return List::create();
}

























