// This file is dependent on fitObj.hpp


// param has 4 rows, a, b, c, lm1.
template <int method, int RIBlib>
NumericVector distances(List distlist, NumericMatrix param, 
                        double eps, int maxCore, int maxit)
{
  int Ndist = distlist.size();
  struct Dist { int size; double *val, *P; };
  if (Ndist != param.ncol()) stop("distlist.size() != param.ncol()");
  std::vector<Dist> dlist(distlist.size());
  for (int i = 0; i < Ndist; ++i)
  {
    List tmp = distlist[i];
    NumericVector val = tmp[0], P = tmp[1];
    dlist[i].val = &val[0];
    dlist[i].P = &P[0];
    dlist[i].size = val.size();
  }
  maxCore = std::min<int> (maxCore, distlist.size());
  CharlieThreadPool cp(maxCore);
  NumericVector result(Ndist);
  std::vector<ObjFun<method, RIBlib> > objfs(maxCore);
  auto f = [&](std::size_t i, std::size_t t)->bool
  {
    objfs[t].reset(dlist[i].val, dlist[i].P, dlist[i].size,
                   param[i * 4 + 3], eps, maxit);
    result[i] = objfs[t].f(&param[i * 4]);
    return false;
  };
  cp.parFor(0, Ndist, f, Ndist / (maxCore * maxCore) + 1);
  return result;
} 




template <int RIBlib>
NumericVector distances(List distlist, NumericMatrix param, 
                        double scaleEps, 
                        int scaleMaxit, 
                        int maxCore,
                        String distanceFun)
{
  NumericVector rst;
  
  
  if (distanceFun == "likelihood") 
    rst = distances<0, RIBlib> (distlist, param, scaleEps, maxCore, scaleMaxit);
  
  
  else if (distanceFun == "likelihoodDiscrete") 
    rst = distances<1, RIBlib> (distlist, param, scaleEps, maxCore, scaleMaxit);
    
  
  else if (distanceFun == "Komogorov") 
    rst = distances<2, RIBlib> (distlist, param, scaleEps, maxCore, scaleMaxit);
    
  
  else if (distanceFun == "EuclideanCDF") 
    rst = distances<3, RIBlib> (distlist, param, scaleEps, maxCore, scaleMaxit);
  
  
  else if (distanceFun == "unbiasedLikelihoodDiscrete") 
    rst = distances<4, RIBlib> (distlist, param, scaleEps, maxCore, scaleMaxit);
  
  
  else if (distanceFun == "unbiasedKomogorov") 
    rst = distances<5, RIBlib> (distlist, param, scaleEps, maxCore, scaleMaxit);
  
  
  else if (distanceFun == "unbiasedEuclideanCDF") 
    rst = distances<6, RIBlib> (distlist, param, scaleEps, maxCore, scaleMaxit);
  
  
  else stop("Distance function not implemented");
  
  
  return rst;
}
  



// [[Rcpp::export]]
NumericVector distances(List distlist, 
                        NumericMatrix param, 
                        double scaleEps = 1e-8, 
                        int scaleMaxit = 100, 
                        int maxCore = 1000,
                        String distanceFun = "likelihood",
                        String RIBlib = "Numerical Recipes")
{
  if (RIBlib == "Boost") 
    return distances<0>(distlist, param, scaleEps, scaleMaxit, maxCore, distanceFun);
  
  
  if (RIBlib == "R::pbeta") 
    return distances<1>(distlist, param, scaleEps, scaleMaxit, maxCore, distanceFun);
  
  
  if (RIBlib == "Numerical Recipes")
    return distances<2>(distlist, param, scaleEps, scaleMaxit, maxCore, distanceFun);
  
  
  stop("Method of computing Regularized Incomplete Beta function is not implemented.");
  return NumericVector(0);
}





































