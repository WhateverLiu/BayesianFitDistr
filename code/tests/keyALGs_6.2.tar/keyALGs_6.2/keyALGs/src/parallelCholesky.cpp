
// #pragma once


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "h/dnyTasking.hpp"
#include <mutex>
using namespace Rcpp;


template<typename num>
bool computeDiag(num *A, std::size_t dim, num *L, std::size_t whichRow)
{
  num *l = L + whichRow * dim;
  std::size_t diagOffset = dim * whichRow + whichRow;
  num tmp = A[diagOffset] - std::inner_product(l, l + whichRow, l, 0.0);
  if(tmp < 0) return false;
  L[diagOffset] = std::sqrt(tmp);
  return true;
}


// Pretend L is row-major. i > j.
template<typename num>
void computeOffDiag(num *A, std::size_t dim, num *L,
                    std::size_t i, std::size_t j)
{
  std::size_t i_dim = i * dim, j_dim = j * dim;
  L[i_dim + j] = (A[j_dim + i] - std::inner_product(
      L + j_dim, L + j_dim + j, L + i_dim, 0.0)) / L[j_dim + j];
}


template<typename num>
struct ComputeOffDiagCol
{
  std::size_t dim, whichCol;
  num *A, *L;
  ComputeOffDiagCol(){}
  ComputeOffDiagCol(std::size_t dim, std::size_t whichCol, num *A, num *L):
    dim(dim), whichCol(whichCol), A(A), L(L){}
  std::size_t Njobs()
  {
    return dim > whichCol + 1 ? dim - whichCol - 1 : 0;
  }
  void run(std::size_t i)
  {
    computeOffDiag(A, dim, L, i + 1 + whichCol, whichCol);
  }
};


/*
template<typename num>
struct ComputeOffDiagCol: public RcppParallel::Worker
{
  std::size_t dim, whichCol;
  num *A, *L;
  dynamicTasking *dT;
  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      std::size_t gap = 1;
      if(!dT->nextTaskID(objI, gap)) break;
      for(std::size_t I = objI, Iend = std::min(objI + gap, dT->NofAtom);
          I < Iend; ++I)
      {
        computeOffDiag(A, dim, L, I + 1 + whichCol, whichCol);
      }
    }
  }
  ComputeOffDiagCol(std::size_t dim, std::size_t whichCol,
                    num *A, num *L, std::size_t maxCore):
    dim(dim), whichCol(whichCol), A(A), L(L)
  {
    if(dim <= whichCol + 1) return;
    dynamicTasking dt(maxCore, dim - whichCol - 1); dT = &dt;
    RcppParallel::parallelFor(0, maxCore, *this);
  }
};
*/




// A is the symmetric matrix to be factorized.
// U is the upper triangle such that U'U = A.
template<typename num>
bool cholesky(num *A, std::size_t dim, num *U, std::size_t maxCore)
{
  if(A[0] < 0) return false;
  U[0] = std::sqrt(A[0]);
  ComputeOffDiagCol<num> cod(dim, 0, A, U);
  std::size_t grainSize = 1;
  ParaFor<ComputeOffDiagCol<num> >(cod, maxCore, grainSize);
  for(std::size_t whichCol = 1; whichCol < dim; ++whichCol)
  {
    bool decomposable = computeDiag(A, dim, U, whichCol);
    if(!decomposable) return false;
    cod.whichCol = whichCol;
    ParaFor<ComputeOffDiagCol<num> >(cod, maxCore, grainSize);
  }
  return true;
}




template<typename num>
struct ComputeOffDiagColFullPara: public RcppParallel::Worker
{
  std::atomic<signed char> *wait, *decomposable;
  std::atomic<std::size_t> *rowCounter, *whichCol;
  std::size_t dim;
  num *A, *L;
  std::size_t maxCore;
  // std::mutex mtx;
  void operator() (std::size_t st, std::size_t end)
  {
    while(*whichCol < dim and *decomposable)
    {
      std::size_t whichRow = rowCounter->fetch_add(1);
      if(st != 0)
      {
        if(whichRow >= dim)
        {
          wait[st] = true;
          while(wait[st]){}
        }
        else computeOffDiag<num> (A, dim, L, whichRow, *whichCol);
      }
      else
      {
        if(whichRow >= dim)
        {
          wait[st] = true;
          while(std::accumulate(wait, wait + maxCore, 0) - maxCore != 0){}
          ++*whichCol;
          if(*whichCol < dim) *decomposable = computeDiag(A, dim, L, *whichCol);
          *rowCounter = *whichCol + 1;
          std::fill(wait, wait + maxCore, false);
        }
        else computeOffDiag<num> (A, dim, L, whichRow, *whichCol);
      }
    }
  }


  ComputeOffDiagColFullPara(std::size_t dim, num *A, num *L, std::size_t maxCore,
                            bool &isDecomposed): dim(dim), A(A),
                            L(L), maxCore(maxCore)
  {
    if(A[0] < 0) { isDecomposed = false; return; }
    L[0] = std::sqrt(A[0]);
    std::atomic<signed char> decomposed(true); decomposable = &decomposed;
    std::atomic<std::size_t> rowCounter_(1); rowCounter = &rowCounter_;
    std::atomic<std::size_t> whichCol_(0); whichCol = &whichCol_;
    std::atomic<signed char> wait_[300]; wait = wait_;
    std::fill(wait, wait + maxCore, false);
    RcppParallel::parallelFor(0, maxCore, *this);
    isDecomposed = decomposed;
  }
};




// [[Rcpp::export]]
NumericMatrix cholDecomp(NumericMatrix A, int maxCore = 15,
                         bool useFullPara = true)
{
  std::size_t dim = A.nrow();
  if(dim - A.ncol() != 0) stop("Not a square matrix.");
  NumericMatrix U(dim, dim);
  if(!useFullPara)
  {
    bool rst = cholesky<double> (&A[0], dim, &U[0], maxCore);
    if(!rst) stop("Matrix is not positive definite.");
    return U;
  }
  bool rst = false;
  ComputeOffDiagColFullPara<double> (dim, &A[0], &U[0], maxCore, rst);
  if(!rst) stop("Matrix is not positive definite.");
  return U;
}



































