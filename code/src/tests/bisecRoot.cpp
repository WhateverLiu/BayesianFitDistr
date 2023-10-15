// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include "../h/bisecRoot.hpp"
using namespace Rcpp;


struct Log
{
  double operator()(double x) { return std::log(x); }
};


// [[Rcpp::export]]
double bisecRootTest(double x, double eps = 1e-6) 
{
  Log f;
  return bisecRoot(x, f, eps);
}
















