// #include <Rcpp.h>
// using namespace Rcpp;



template <int RIBlib>
NumericVector solve_d(NumericVector abc_lm1, double eps, int maxit, 
                      bool useNewton) 
{
  TrB trb;
  double *param = &abc_lm1[0];
  NumericVector rst(abc_lm1.size() / 4);
  for (int i = 0, iend = abc_lm1.size(); i < iend; i += 4)
  {
    double *p = param + i;
    if (useNewton) trb.reset<true, RIBlib> (p[0], p[1], p[2], p[3], eps, maxit);
    else trb.reset<false, RIBlib> (p[0], p[1], p[2], p[3], eps, maxit);
    rst[i >> 2] = trb.d;
  }
  return rst;
}




// [[Rcpp::export]]
NumericVector solve_d(NumericVector abc_lm1, double eps, int maxit, 
                      bool useNewton = true, String RIBlib = "R::pbeta")
{
  if (RIBlib == "Boost") return solve_d<0> (abc_lm1, eps, maxit, useNewton);
  if (RIBlib == "R::pbeta") return solve_d<1> (abc_lm1, eps, maxit, useNewton);
  if (RIBlib == "Numerical Recipes") return solve_d<2> (abc_lm1, eps, maxit, useNewton);
  stop("Method of computing Regularized Incomplete Beta function is not implemented.");
  return NumericVector(0);
}
  




