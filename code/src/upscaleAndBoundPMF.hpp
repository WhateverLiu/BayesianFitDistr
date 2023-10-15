template <int regridMethod>
List upScaleAndBoundPMF(List pmf, double upperBound, NumericVector upscaler,
                        bool regridedMaxUnchanged = true)
{
  List rst(upscaler.size());
  NumericVector val = pmf[0], P = pmf[1];
  
  
  // for (int k = 1, kend = val.size(); k < kend; ++k)
  //   if (val[k] <= val[k - 1]) stop("val[k] <= val[k - 1]!");
  
  
  // int sizeold = val.size();
  std::vector<double> cntr(val.size() * 2);
  Regrid regrider;
  
  
  for (int i = 0, iend = upscaler.size(); i < iend; ++i)
  {
    double r = upscaler[i];
    int sizenew = val.size();
    bool lastValShouldBeUpperBound = false;
    for (int j = 0, jend = val.size(); j < jend; ++j)
    {
      double v = val[j] * r;
      if (v >= upperBound) { lastValShouldBeUpperBound = true; sizenew = j + 1; break; }
    }
    
    
    double *valnew = cntr.data(), *Pnew = valnew + sizenew;
    double psum = 0;
    for (int j = 0, jend = sizenew - 1; j < jend; ++j)
    {
      valnew[j] = val[j] * r;
      Pnew[j] = P[j];
      psum += P[j];
    }
    if (lastValShouldBeUpperBound) valnew[sizenew - 1] = upperBound;
    else valnew[sizenew - 1] = val[sizenew - 1] * r;
    Pnew[sizenew - 1] = 1.0 - psum;
    
    
    NumericVector u(sizenew), v(sizenew);
    double valmax = valnew[sizenew - 1];
    if (r < 1 and regridedMaxUnchanged) valmax = val[val.size() - 1];
    makeRegularGrid(valnew[0], valmax, &u[0], sizenew);
    regrider.operator()<
      regridMethod, int, double, double, double, double, false, true> (
          valnew, Pnew, sizenew, &u[0], &v[0], sizenew);
    
    
    rst[i] = List::create(Named("val") = u, Named("P") = v);
  }
  return rst;
}


// [[Rcpp::export]]
List upScaleAndBoundPMF(List pmf, double upperBound, NumericVector upscaler,
                        String regridMethod = "lmm")
{
  if (regridMethod == "lr")
    return upScaleAndBoundPMF<0>(pmf, upperBound, upscaler);
  if (regridMethod == "lmm")
    return upScaleAndBoundPMF<1>(pmf, upperBound, upscaler);
  if (regridMethod == "r4")
    return upScaleAndBoundPMF<2>(pmf, upperBound, upscaler);
  stop("Regrid method not implemented.");
  return List::create();
}











