// #pragma once
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
# include <Rcpp.h>
// # include <RcppArmadillo.h>
# include <RcppParallel.h>
# include "h/dnyTasking.hpp"
// # include <atomic>
using namespace Rcpp;
// # define vec std::vector


template<typename num>
inline void vecxscalerAdd(num *x, std::size_t size, num scalar, num *rst)
{
  for(std::size_t i = 0; i < size; ++i)
    rst[i] += x[i] * scalar;
}


template<typename num>
inline void vecxscalerAddBackward(num *x, std::size_t size, num scalar, num *rst)
{
  for(std::size_t i = size - 1, iend = 0 - 1; i != iend; --i)
    rst[i] += x[i] * scalar;
}


// y is a column vector and of size ncol.
template<typename num>
inline void matxcol(num *X, std::size_t nrow, std::size_t ncol, num *y, num *rst)
{
  std::fill(rst, rst + nrow, 0);
  for(std::size_t i = 0; i < ncol; ++i)
    vecxscalerAdd(X + i * nrow, nrow, y[i], rst);
}


template<typename num>
inline void matxcolBackward(num *X, std::size_t nrow, std::size_t ncol, num *y, num *rst)
{
  std::fill(rst, rst + nrow, 0);
  for(std::size_t i = ncol - 1, iend = 0 - 1; i != iend; --i)
    vecxscalerAddBackward(X + i * nrow, nrow, y[i], rst);
}


template<typename num>
inline void matmulSingleThreadBidir(num *X, std::size_t nrowX, std::size_t ncolX,
                               num *Y, std::size_t ncolY, num *rst)
{
  for(std::size_t i = 0; i < ncolY; ++i)
  {
    matxcol(X, nrowX, ncolX, Y + ncolX * i, rst + nrowX * i);
    ++i;
    if(i >= ncolY) break;
    matxcolBackward(X, nrowX, ncolX, Y + ncolX * i, rst + nrowX * i);
  }
}


template<typename num>
inline void matmulSingleThread(num *X, std::size_t nrowX, std::size_t ncolX,
                               num *Y, std::size_t ncolY, num *rst)
{
  for(std::size_t i = 0; i < ncolY; ++i)
  {
    matxcol(X, nrowX, ncolX, Y + ncolX * i, rst + nrowX * i);
  }
}


template<typename ing, typename num>
struct paraMatMul: public RcppParallel::Worker
{
  ing grainSize;
  ing N, P; // N x P, P x M
  num *X, *Y, *rst;
  dynamicTasking *dT;
  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI, grainSize)) break; // Number of tasks equals M
      for(std::size_t i = objI, iend = std::min(i + grainSize, dT->NofAtom);
          i < iend; ++i)
      {
        matxcol(X, N, P, Y + P * i, rst + N * i);
      }
    }
  }
  paraMatMul(num *X, num *Y, num *rst,
             ing N, ing P, ing M, ing maxCore): // (N x P) x (P x M)
    N(N), P(P), X(X), Y(Y), rst(rst)
  {
    maxCore = std::min<ing> (maxCore, M);
    dynamicTasking dt(maxCore, M); dT = &dt;
    grainSize = std::max<ing> (1, M / (maxCore * maxCore));
    parallelFor(0, maxCore, *this);
  }
};








// [[Rcpp::export]]
NumericMatrix matmul(NumericMatrix X, NumericMatrix Y, int maxCore = 16,
                     bool bidir = false)
{
  NumericMatrix rst(X.nrow(), Y.ncol());
  if(maxCore == 1)
  {
    if(!bidir) matmulSingleThread<double> (
        &X[0], X.nrow(), X.ncol(), &Y[0], Y.ncol(), &rst[0]);
    else matmulSingleThreadBidir<double> (
        &X[0], X.nrow(), X.ncol(), &Y[0], Y.ncol(), &rst[0]);
  }
  else
    paraMatMul<int, double> (
        &X[0], &Y[0], &rst[0], X.nrow(), X.ncol(), Y.ncol(), maxCore);
  return rst;
}


// Faster than armadillo actually.
// // [[Rcpp::export]]
// arma::mat matmulArmadillo(arma::mat &X, arma::mat &Y)
// {
//   return X * Y;
// }























