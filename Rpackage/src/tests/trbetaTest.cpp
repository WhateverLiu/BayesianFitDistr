// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include "../trb/trb.hpp"
using namespace Rcpp;


// [[Rcpp::export]]
List testTrbLimMean(List x, List abc)
{
  List rst(x.size());
  TrbLimMean t;
  for (int i = 0, iend = rst.size(); i < iend; ++i)
  {
    NumericVector d = x[i];
    NumericVector p = abc[i];
    t.reset(p[0], p[1], p[2]);
    NumericVector rsti(d.size());
    for (int j = 0, jend = d.size(); j < jend; ++j)
      rsti[j] = t(d[j]);
    rst[i] = rsti;
  }
  return rst;
}


// [[Rcpp::export]]
List testTrbLogPdf(List x, List abcd)
{
  List rst(x.size());
  TrbLogPdf t;
  for (int i = 0, iend = rst.size(); i < iend; ++i)
  {
    NumericVector d = x[i];
    NumericVector p = abcd[i];
    t.reset(p[0], p[1], p[2], p[3]);
    NumericVector rsti(d.size());
    for (int j = 0, jend = d.size(); j < jend; ++j)
      rsti[j] = t(d[j]);
    rst[i] = rsti;
  }
  return rst;
}


// [[Rcpp::export]]
List testTrbCdf(List x, List abcd)
{
  List rst(x.size());
  TrbCdf t;
  for (int i = 0, iend = rst.size(); i < iend; ++i)
  {
    NumericVector d = x[i];
    NumericVector p = abcd[i];
    t.reset(p[0], p[1], p[2], p[3]);
    NumericVector rsti(d.size());
    for (int j = 0, jend = d.size(); j < jend; ++j)
      rsti[j] = t(d[j]);
    rst[i] = rsti;
  }
  return rst;
}















// [[Rcpp::export]]
double testlevtrbeta(
    double shape1, double shape2, double shape3,
    double scale)
{
  return levtrbeta(shape1, shape2, shape3, scale);
}


// [[Rcpp::export]]
double testDtrbeta(
    double x, double shape1, double shape2, double shape3,
    double scale, bool give_log)
{
  return dtrbeta(x, shape1, shape2, shape3,
                 scale, give_log);
}


// [[Rcpp::export]]
double testPtrbeta(
    double q, double shape1, double shape2, double shape3,
    double scale, int lower_tail, int log_p)
{
  return ptrbeta(q, shape1, shape2, shape3,
                 scale, lower_tail, log_p);
}












































































