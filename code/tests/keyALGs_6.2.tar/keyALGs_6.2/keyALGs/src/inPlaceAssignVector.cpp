// #pragma once

#include <Rcpp.h>
using namespace Rcpp;


template<typename T>
void vecAssign(T *x, int xsize, int *ind, T *y, int size)
{
  for(int i = 0; i < size; ++i)
  {
    if(ind[i] >= xsize) stop("Index out of bound.");
    x[ind[i] - 1] = y[i];
  }
}


// [[Rcpp::export]]
void testAddr(IntegerVector x)
{
  Rcout << &x[0] << "\n";
}


// [[Rcpp::export]]
void vecAssign(SEXP orginalV, IntegerVector index, SEXP newValueV)
{
  if(TYPEOF(orginalV) == 13)
  {
    IntegerVector x = orginalV, y = newValueV;
    vecAssign<int> (&x[0], x.size(), &index[0], &y[0], index.size());
  }
  else if(TYPEOF(orginalV) == 14)
  {
    NumericVector x = orginalV, y = newValueV;
    vecAssign<double> (&x[0], x.size(), &index[0], &y[0], index.size());
  }
  else if(TYPEOF(orginalV) == 24)
  {
    RawVector x = orginalV, y = newValueV;
    vecAssign<unsigned char> (&x[0], x.size(), &index[0], &y[0], index.size());
  }
  else stop("value type not implemented.");
}


/***R
if(F)
{
  tmp = as.integer(runif(100000) * 1e6)
  testAddr(tmp)
  vecAssign(tmp, c(3, 2), c(10000, 20000))
  testAddr(tmp)
}
*/


