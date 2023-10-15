// [[Rcpp::plugins(cpp17)]]
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// #include <Rcpp.h>
using namespace Rcpp;
#include "../lbfgs/LBFGSBcharlie.hpp"
#include "../h/plainGradientDescent.hpp"



struct Rosenbrock
{
  double operator()(double *x, double *&grad, int dim)
  {
    // grad.resize(dim);
    double fx = 0.0;
    for(int i = 0; i < dim; i += 2)
    {
      double t1 = 1.0 - x[i];
      double t2 = 10 * (x[i + 1] - x[i] * x[i]);
      grad[i + 1] = 20 * t2;
      grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
      fx += t1 * t1 + t2 * t2;
    } 
    return fx;
  }
};


struct RosenbrockNoGrad
{
  double operator()(double *x, double *&grad, int dim)
  {
    // grad.resize(dim);
    double fx = 0.0;
    for(int i = 0; i < dim; i += 2)
    { 
      double t1 = 1.0 - x[i];
      double t2 = 10 * (x[i + 1] - x[i] * x[i]);
      fx += t1 * t1 + t2 * t2;
    }
    grad = nullptr;
    return fx;
  } 
};






// [[Rcpp::export]]
List miniRosenbrock(NumericVector x, int maxit = 100, bool giveGrad = false,
                    NumericVector lb = NumericVector(0), 
                    NumericVector ub = NumericVector(0)) 
{
  if (x.size() % 2 != 0) stop("x.size() is not even!");
  NumericVector rst(x.begin(), x.end());
  LBFGSBcharlie lbfgsb;
  Rcout << "sizeof(lbfgsb) = " << sizeof(lbfgsb) << "\n";
  
  
  double *lbptr = nullptr, *ubptr = nullptr;
  if (lb.size() != 0) lbptr = &lb[0];
  if (ub.size() != 0) ubptr = &ub[0];
  
  
  if (giveGrad)
  {
    Rosenbrock f;
    auto fval = lbfgsb(&rst[0], x.size(), f, maxit, lbptr, ubptr);
    return List::create(Named("param") = rst, Named("fval") = fval.first,
                        Named("Niter") = fval.second);  
  }
  
  
  RosenbrockNoGrad f;
  auto fval = lbfgsb(&rst[0], x.size(), f, maxit, lbptr, ubptr);
  return List::create(Named("param") = rst, Named("fval") = fval.first,
                      Named("Niter") = fval.second);  
  
}




// [[Rcpp::export]]
List miniRosenbrockPlainGD(
    NumericVector x, int maxit = 100, bool giveGrad = false, double lr = 1e-4) 
{
  if (x.size() % 2 != 0) stop("x.size() is not even!");
  NumericVector rst(x.begin(), x.end());
  
  
  PlainGD gd;
  if (giveGrad)
  {
    Rosenbrock f;
    auto fval = gd(&rst[0], x.size(), f, lr, maxit);
    return List::create(Named("param") = rst, Named("fval") = fval.first,
                        Named("Niter") = fval.second);  
  }
  
  
  RosenbrockNoGrad f;
  auto fval = gd(&rst[0], x.size(), f, lr, maxit);
  return List::create(Named("param") = rst, Named("fval") = fval.first,
                      Named("Niter") = fval.second);  
  
}












