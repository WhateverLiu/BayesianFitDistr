// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include "../trb/trb.hpp"
using namespace Rcpp;


// [[Rcpp::export]]
double levtrbetaTest ( double shape1, double shape2, double shape3,
                 double scale) {
  // return R::pbeta(x1m, bp, ap, /*l._t.*/0, /*give_log*/1) 
  // + Rmath::gammafun(x1m)
  // ;
  return levtrbeta(// double limit, 
    shape1, shape2, shape3,
    scale //, double order, int give_log
  );
}
















