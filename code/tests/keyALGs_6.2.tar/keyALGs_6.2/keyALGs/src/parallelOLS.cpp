// #pragma once


// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>
# include <RcppParallel.h>
// # define ARMA_USE_SUPERLU 1
// # define ARMA_SUPERLU_INCLUDE_DIR SuperLU/SRC/
# include <time.h>
using namespace Rcpp;
using namespace RcppParallel;
# define valtype double
# define indtype int
# define vec std::vector
# define INT std::size_t
# define eps 1e-32
# define Inf 1e308


struct dynamicTasking
{
  std::size_t NofCore;
  std::size_t NofAtom;
  tbb::atomic<std::size_t> counter;


  void reset(std::size_t NofCPU, std::size_t NofTask)
  {
    NofCore = NofCPU;
    if(NofCore > NofTask) NofCore = NofTask;
    NofAtom = NofTask;
    counter = 0;
  }


  dynamicTasking(){}


  dynamicTasking(std::size_t NofCPU, std::size_t NofTask)
  {
    reset(NofCPU, NofTask);
  }


  bool nextTaskID(std::size_t &taskID)
  {
    taskID = counter.fetch_and_increment();
    return taskID < NofAtom;
  }


  bool nextTaskID(std::size_t &taskID, std::size_t increment)
  {
    taskID = counter.fetch_and_add(increment);
    return taskID < NofAtom;
  }
};




/*
// ============================================================================
// BEGIN:
// Parallel Gaussian elimination to triangle matrix U then back-substitution
// ============================================================================


// X is the data matrix of size d x d, y is a vector size d.
// X stores the data column by column.
// rst is resulting vector of size d.
template<bool numericHealthy>
struct oneGaussianEl: public Worker
{
  indtype dim; // dimension
  indtype whichRowIsAddedToRowsBelow;
  valtype **cols;
  valtype *multipliers;
  dynamicTasking *dT;
  void operator() (INT st, INT end)
  {
    for(;;)
    {
      INT objI = 0;
      if(!dT->nextTaskID(objI, 256)) break;
      for(INT I = objI, Iend = std::min<INT> (I + 256, dT->NofAtom);
          I < Iend; ++I)
      {
        valtype *v = cols[whichRowIsAddedToRowsBelow + I + 1];
        for(indtype i = whichRowIsAddedToRowsBelow + 1; i < dim; ++i)
        {
          if(numericHealthy) v[i] += v[whichRowIsAddedToRowsBelow] * multipliers[i];
          // else
          // {
          //   valtype tmp = v[i] + v[whichRowIsAddedToRowsBelow] * multipliers[i];
          //   int sn = std::signbit(tmp) * int(-2) + 1;
          //   tmp = std::min<double> (std::abs(tmp), Inf) * sn;
          //   v[i] = tmp;
          // }
        }
      }
    }
  }


  oneGaussianEl(indtype dim, indtype whichRowIsAddedToRowsBelow,
                valtype **cols,
                valtype *multipliers, // A container of size dim
                int maxCore):
    dim(dim), whichRowIsAddedToRowsBelow(whichRowIsAddedToRowsBelow), cols(cols),
    multipliers(multipliers)
  {
    valtype *v = cols[whichRowIsAddedToRowsBelow];
    for(indtype i = whichRowIsAddedToRowsBelow + 1; i < dim; ++i)
    {
      if(numericHealthy) multipliers[i] = -v[i] / v[whichRowIsAddedToRowsBelow];
      // else
      // {
      //   valtype tmp = -v[i] / v[whichRowIsAddedToRowsBelow];
      //   int sn = std::signbit(tmp) * (-2) + 1;
      //   tmp = std::min<double> (std::abs(tmp), Inf) * sn;
      //   multipliers[i] = tmp;
      // }
    }
    // std::cout << "In paralle for constructor\n";
    dT = new dynamicTasking(maxCore, dim - whichRowIsAddedToRowsBelow);
    parallelFor(0, maxCore, *this);
    delete dT;
  }
};


// Contents in cols will be contaminated
void gaussianToU(valtype **cols, indtype dim, // dim is the number of rows
                 bool numericHealty, int maxCore)
{
  vec<valtype> multipliers(dim);
  // auto start = std::chrono::system_clock::now();
  for(indtype i = 0, iend = dim - 1; i < iend; ++i)
  {
    indtype &whichRowIsAddedToRowsBelow = i;
    if(numericHealty) oneGaussianEl<true> (
          dim, whichRowIsAddedToRowsBelow, cols, &multipliers[0], maxCore);
    else oneGaussianEl<false> (
        dim, whichRowIsAddedToRowsBelow, cols, &multipliers[0], maxCore);
  }
  // auto tend = std::chrono::system_clock::now();
  // std::chrono::duration<double> diff = tend - start;
  // std::cout << "timeSpent = " << diff.count() << "\n";
}


// [[Rcpp::export]]
NumericMatrix gaussianEliminationToUpperTriangle(
    NumericMatrix X, NumericVector y, bool numericHealty = true, int maxCore = 7)
{
  indtype dim = X.ncol();
  NumericMatrix rst(dim, dim + 1);
  std::copy(X.begin(), X.end(), rst.begin());
  std::copy(y.begin(), y.end(), rst.end() - dim);
  vec<valtype*> cols(dim + 1);
  for(indtype i = 0; i < dim + 1; ++i)
  {
    cols[i] = &rst[0] + INT(dim) * i;
  }
  gaussianToU(&cols[0], dim, numericHealty, maxCore);
  return rst;
}



template<bool numericHealthy>
struct backSubstituteOneRow: public Worker
{
  indtype dim;
  indtype solveWhichRow;
  valtype *X; // of dim rows and dim columns
  valtype *rst;
  valtype *threadSum;
  dynamicTasking *dT;
  void operator() (INT st, INT end)
  {
    for(;;)
    {
      INT objI = 0;
      if(!dT->nextTaskID(objI, 64)) break;
      for(INT I = objI, // i is target's column index
          Iend = std::min<INT> (I + 64, dT->NofAtom);
          I < Iend; ++I)
      {
        INT i = I + solveWhichRow + 1;
        threadSum[st] += rst[i] * X[dim * i + solveWhichRow];
      }
    }
  }


  backSubstituteOneRow(bool &unsolvable, indtype dim, indtype solveWhichRow,
                       valtype *X, valtype *y, valtype *rst,
                       vec<valtype> &threadSumContainer, int maxCore):
    dim(dim), solveWhichRow(solveWhichRow), X(X), rst(rst)
  {
    threadSumContainer.assign(maxCore, 0);
    threadSum = &threadSumContainer[0];
    dT = new dynamicTasking(maxCore, dim - solveWhichRow - 1);
    parallelFor(0, maxCore, *this);
    delete dT;
    valtype S = std::accumulate(threadSum, threadSum + maxCore, 0.0);
    valtype coef = X[INT(dim) * solveWhichRow + solveWhichRow];
    if(coef == 0)
    {
      unsolvable = true;
      return;
    }


    if(numericHealthy) rst[solveWhichRow] = (y[solveWhichRow] - S) / coef;
    else
    {
      // int sn = std::signbit(S) * (-2) + 1;
      // S = std::min<double> (std::abs(S), Inf) * sn;
    }
  }
};


// X has been trianglized
bool paraSolve(valtype *rst, valtype *X, valtype *y, indtype dim,
               bool numericHealty, int maxCore)
{
  bool unsolvable = false;
  vec<valtype*> cols(dim + 1);
  for(indtype i = 0, iend = dim + 1; i < iend; ++i)
    cols[i] = X + INT(dim) * i;
  cols.back() = y;
  gaussianToU(&cols[0], dim, numericHealty, maxCore);


  vec<valtype> threadSumContainer(maxCore);
  for(indtype i = dim - 1; i >= 0; --i)
  {
    indtype &solveWhichRow = i;
    if(numericHealty)
      backSubstituteOneRow<true> (
          unsolvable, dim, solveWhichRow, X, y, rst, threadSumContainer, maxCore);
    else
      backSubstituteOneRow<false> (
          unsolvable, dim, solveWhichRow, X, y, rst, threadSumContainer, maxCore);
    if(unsolvable) return true;
  }
  return false;
}


// [[Rcpp::export]]
arma::colvec paraSolve(NumericMatrix X, arma::colvec y, bool numericHealty = true,
                       int maxCore = 7)
{
  indtype dim = X.ncol();
  arma::colvec rst(dim);
  arma::mat Xmat(dim, dim);
  std::copy(X.begin(), X.end(), Xmat.begin());
  bool unsolvable = paraSolve(
    &*rst.begin(), &*Xmat.begin(), &*y.begin(), dim, numericHealty, maxCore);
  if(unsolvable) return arma::colvec(0.0);
  return rst;
}





// ============================================================================
// END:
// Parallel Gaussian elimination to triangle matrix U then back-substitution
// ============================================================================
*/




// [[Rcpp::export]]
arma::colvec armaSolve(arma::mat &X, arma::colvec &y, bool fastsolve = true)
{
  if(fastsolve) return arma::solve(X, y, arma::solve_opts::fast);
  else return arma::solve(X, y);
}




struct couple
{
  indtype j, i;
};
struct XtX: public Worker
{
  const unsigned verboseInterval;
  indtype nrow, ncol;
  couple *pairInd;
  valtype **col;
  valtype **rst; // rst points to a container of ncol x ncol
  dynamicTasking *dT;
  void operator()(std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI)) break;
      {
        if((objI & (verboseInterval - 1)) == 0L)
        {
          std::cout << 0 - (int64_t)objI;
        }
        indtype &i = pairInd[objI].i, &j = pairInd[objI].j;
        rst[j][i] = std::inner_product(col[i], col[i] + nrow, col[j], 0.0);
        rst[i][j] = rst[j][i];
      }
    }
  }
  XtX(unsigned verboseInterval, indtype nrow, indtype ncol,
      valtype *X, valtype *result, int maxCore):
    verboseInterval(verboseInterval), nrow(nrow), ncol(ncol)
  {
    vec<valtype*> Xcntr(ncol);
    for(indtype i = 0; i < ncol; ++i) Xcntr[i] = X + nrow * i;


    vec<valtype*> resultCntr(ncol);
    for(indtype i = 0; i < ncol; ++i) resultCntr[i] = result + ncol * i;


    col = &Xcntr[0];
    rst = &resultCntr[0];
    vec<couple> pairs(std::size_t(ncol) * (ncol + 1) / 2);
    std::size_t k = 0;
    for(indtype i = 0; i < ncol; ++i)
    {
      for(indtype j = i; j < ncol; ++j)
      {
        pairs[k].j = j;
        pairs[k].i = i;
        ++k;
      }
    }
    pairInd = &pairs[0];
    dynamicTasking dt(maxCore, pairs.size());
    dT = &dt;
    Rcout << "Compute X'X, " << pairs.size() << " inner products:\n";
    parallelFor(0, maxCore, *this);
    Rcout << "\n";
  }
};




struct XtY: public Worker
{
  indtype nrow, ncol;
  valtype *x, *y, *xty;
  dynamicTasking *dT;
  void operator()(std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI)) break;
      {
        xty[objI] = std::inner_product(
          x + objI * nrow, x + (objI + 1) * nrow, y, 0.0);
      }
    }
  }
  XtY(indtype nrow, indtype ncol,
      valtype *x, valtype *y, valtype *xty, int maxCore):
    nrow(nrow), ncol(ncol), x(x), y(y), xty(xty)
  {
    dynamicTasking dt(maxCore, ncol); dT = &dt;
    parallelFor(0, maxCore, *this);
  }
};


struct calPred: public Worker
{
  indtype nrow, ncol;
  valtype *x, *w;
  vec<valtype> *tmp;
  dynamicTasking *dT;
  void operator()(std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI)) break;
      {
        valtype *v = x + objI * nrow;
        for(indtype i = 0; i < nrow; ++i)
        {
          tmp[st][i] += v[i] * w[objI];
        }
      }
    }
  }
  calPred(indtype nrow, indtype ncol, valtype *x, valtype *w, valtype *rst, int maxCore):
    nrow(nrow), ncol(ncol), x(x), w(w)
  {
    vec<vec<valtype> > tmpCntr(maxCore, vec<valtype>(nrow, 0));
    tmp = &tmpCntr[0];
    dynamicTasking dt(maxCore, ncol);
    dT = &dt;
    parallelFor(0, maxCore, *this);
    for(indtype i = 0; i < maxCore; ++i)
    {
      for(indtype j = 0; j < nrow; ++j)
      {
        rst[j] += tmp[i][j];
      }
    }
  }
};


valtype relaErr(valtype preval, valtype val)
{
  if(preval == 0) return val - preval;
  return val / preval - 1;
}


// [[Rcpp::export]]
arma::mat guassSeidelSolveDense(
    arma::mat &X, arma::colvec &y,
    int maxIter = 1e6, double convergenceEPS = 1e-6)
{
  indtype dim = y.size();
  vec<valtype> ycopy(y.begin(), y.end());
  vec<valtype> Xcopy(X.begin(), X.end());
  valtype *x = &Xcopy[0];
  vec<unsigned char> mapped(dim, false);
  vec<indtype> pivots(dim);
  for(indtype i = 0; i < dim; ++i)
  {
    for(indtype j = 0; j < dim; ++j)
    {
      if(mapped[j] or x[j] == 0) continue;
      mapped[j] = true;
      valtype r = -1.0 / x[j];
      for(indtype k = 0; k < dim; ++k)
      {
        x[k] *= r;
      }
      x[j] = 0;
      pivots[i] = j;
      ycopy[i] *= r;
      break;
    }
    x += dim;
  }


  indtype iter = 0;
  valtype maxErr = 0;
  arma::mat solutionContainer(ycopy.size(), 1);
  double *solution = &solutionContainer[0];
  std::copy(ycopy.begin(), ycopy.end(), solution);
  while(iter < maxIter)
  {
    x = &Xcopy[0];
    maxErr = 0;
    for(indtype i = 0; i < dim; ++i, x += dim)
    {
      valtype pre = solution[pivots[i]];
      solution[pivots[i]] = std::inner_product(
        x, x + dim, solution, 0.0) - ycopy[pivots[i]];
      maxErr = std::max<valtype> (
        maxErr, std::abs(relaErr(pre, solution[pivots[i]])));
    }
    if(maxErr < convergenceEPS) break;
    ++iter;
  }
  Rcout << "Iteration = " << iter;
  if(maxErr < convergenceEPS) Rcout << ", converged.\n";
  else Rcout << ", unconverged! Maximal error = " << maxErr << ".\n";
  return solutionContainer;
}


// [[Rcpp::export]]
List paraOLS(NumericVector y, NumericMatrix X,
             int maxCore = 7, int verboseInt = 8192,
             bool useGaussSeidel = false, int maxIter = 1e6,
             double convergenceEPS = 1e-6)
{
  int ncol = X.ncol(), nrow = X.nrow();
  arma::mat XtXmat(ncol, ncol);
  XtX(verboseInt, nrow, ncol, &*X.begin(), &*XtXmat.begin(), maxCore);
  arma::colvec XtYvec(ncol);
  XtY(nrow, ncol, &*X.begin(), &*y.begin(), &*XtYvec.begin(), maxCore);
  if(verboseInt > 0) Rcout << "Solving symmetric linear system...\n";
  arma::mat coef;
  if(!useGaussSeidel) arma::solve(coef, XtXmat, XtYvec);
  else
  {
    arma::mat tmp = guassSeidelSolveDense(XtXmat, XtYvec, maxIter, convergenceEPS);
    coef.swap(tmp);
  }
  NumericVector pred(nrow);
  calPred(nrow, ncol, &*X.begin(), &*coef.begin(), &*pred.begin(), maxCore);
  double totalErr = 0;
  for(indtype i = 0; i < nrow; ++i)
  {
    totalErr += (pred[i] - y[i]) * (pred[i] - y[i]);
  }
  return List::create(Named("coef") = coef, Named("pred") = pred,
                      Named("totalErr") = totalErr);
}








struct sparseV
{
  indtype size;
  indtype *ind;
  valtype *val;
};


valtype innerProd(sparseV &x, sparseV &y)
{
  valtype S = 0;
  indtype i = 0, j = 0;
  while(i < x.size and j < y.size)
  {
    if(x.ind[i] == y.ind[j])
    {
      S += x.val[i] * y.val[j];
      ++i; ++j;
    }
    else if(x.ind[i] > y.ind[j]) ++j;
    else ++i;
  }
  return S;
}


struct XtXsparse: public Worker
{
  const unsigned verboseInterval;
  indtype nrow, ncol;
  couple *pairInd;
  sparseV *col;
  valtype **rst; // rst points to a container of ncol x ncol
  dynamicTasking *dT;
  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI, 64)) break;
      for(std::size_t I = objI, Iend = std::min<std::size_t> (dT->NofAtom, I + 64);
          I < Iend; ++I)
      {
        if((I & (verboseInterval - 1)) == 0L)
        // If verboseInterval is power of 2, I & (verboseInterval - 1) == I % verboseInterval
        {
          std::cout << 0 - (int64_t)I;
        }
        indtype &i = pairInd[I].i, &j = pairInd[I].j;
        rst[j][i] = innerProd(col[i], col[j]);
        rst[i][j] = rst[j][i];
      }
    }
  }


  XtXsparse(unsigned verboseInterval, indtype ncol,
      sparseV *X, valtype *result, int maxCore):
    verboseInterval(verboseInterval), ncol(ncol), col(X)
  {
    vec<valtype*> resultCntr(ncol);
    for(std::size_t i = 0, iend = ncol; i < iend; ++i)
      resultCntr[i] = result + ncol * i;


    rst = &resultCntr[0];
    vec<couple> pairs(std::size_t(ncol) * (ncol + 1) / 2);
    std::size_t k = 0;
    for(indtype i = 0; i < ncol; ++i)
    {
      for(indtype j = i; j < ncol; ++j)
      {
        pairs[k].j = j;
        pairs[k].i = i;
        ++k;
      }
    }


    pairInd = &pairs[0];
    dynamicTasking dt(maxCore, pairs.size());
    dT = &dt;
    Rcout << "Compute X'X, " << pairs.size() << " inner products:\n";
    parallelFor(0, maxCore, *this);
    Rcout << "\n";
  }
};


// [[Rcpp::export]]
List fullMatToSparseCol(NumericMatrix Xmat)
{
  INT ncol = Xmat.ncol(), nrow = Xmat.nrow();
  INT nonZero = 0;
  double *X = &Xmat[0];
  for(INT i = 0, iend = ncol * nrow; i < iend; ++i)
  {
    if(X[i] != 0) ++nonZero;
  }
  IntegerVector indexVec(nonZero);
  int *iv = &indexVec[0];
  NumericVector valVec(nonZero);
  double *vv = &valVec[0];
  IntegerVector siz(ncol);


  for(INT i = 0; i < ncol; ++i)
  {
    valtype *x = X + i * nrow;
    INT k = 0;
    for(INT j = 0; j < nrow; ++j)
    {
      if(x[j] != 0)
      {
        iv[k] = j;
        vv[k] = x[j];
        ++k;
      }
    }
    siz[i] = k;
    iv += k;
    vv += k;
  }
  return List::create(Named("size") = siz, Named("index") = indexVec, Named("value") = valVec);
}




template<typename T>
inline T *adrs(void *current)
{
  return (T*)((INT(current) + (sizeof(INT) - 1)) & ~(sizeof(INT) - 1));
}




/*
// ============================================================================
// BEGIN:
// Parallel Gaussian elimination to triangle matrix U then back-substitution
// for sparse matrix. Tests show parallelizing Gaussian elimination mostly
// backfires because jobs on each thread are extremely simple. Setting maxCore
// to 1 is highly recommended
// ============================================================================
struct sparseGauEli: public Worker
{
  indtype dim;
  indtype whichRowAddsToRowsBelow;
  sparseV *rows;
  valtype *y;
  valtype *solution;
  dynamicTasking *dT;


  void operator() (INT st, INT end)
  {
    for(;;)
    {
      INT objI = 0;
      if(!dT->nextTaskID(objI, 512)) break;
      sparseV *add = rows + whichRowAddsToRowsBelow;
      valtype *addY = y + whichRowAddsToRowsBelow;


      for(INT I = objI, Iend = std::min<INT> (I + 512, dT->NofAtom);
          I < Iend; ++I)
      {
        sparseV *receiver = add + 1 + I;
        valtype *receiverY = addY + 1 + I;
        // std::cout << "*receiverY = " << *receiverY << "\n";
        if(receiver->size == 0) continue;
        if(receiver->size == 1)
        {
          solution[receiver->ind[0]] = *receiverY / receiver->val[0];
          receiver->size = 0;
          continue;
        }


        indtype pivot = add->ind[0];
        // std::cout << "pivot = " << pivot << "\n";
        if(receiver->ind[0] > pivot) continue;
        indtype which = std::lower_bound(receiver->ind, receiver->ind +
          receiver->size, pivot) - receiver->ind;
        if(which >= receiver->size) continue;


        valtype r = -receiver->val[which] / add->val[0];
        // std::cout << "r = " << r << "\n";
        // std::cout << "which = " << which << "\n";
        *receiverY += *addY * r;
        {
          indtype i = 0, j = 0;
          while(i < add->size and j < receiver->size)
          {
            if(add->ind[i] == receiver->ind[j])
            {
              receiver->val[j] += r * add->val[i];
              ++i; ++j;
            }
            else if(add->ind[i] > receiver->ind[j]) ++j;
            else ++i;
          }
        }


        if(which < receiver->size / 2)
        {
          std::copy_backward(receiver->ind, receiver->ind + which, receiver->ind + 1);
          ++receiver->ind;
          std::copy_backward(receiver->val, receiver->val + which, receiver->val + 1);
          ++receiver->val;
        }
        else
        {
          std::copy(receiver->ind + which + 1, receiver->ind + receiver->size,
                    receiver->ind + which);
          std::copy(receiver->val + which + 1, receiver->val + receiver->size,
                    receiver->val + which);
        }
        --receiver->size;
      }
    }
  }


  sparseGauEli(indtype dim, sparseV *rows, valtype *y,
               valtype *solution, int maxCore):
    dim(dim), rows(rows), y(y), solution(solution)
  {
    dynamicTasking dt; dT = &dt;
    for(indtype i = 0, iend = dim - 1; i < 2; ++i)
    {
      if(rows[i].size == 0) continue;
      if(rows[i].size == 1)
      {
        solution[rows[i].ind[0]] = y[i] / rows[i].val[0];
        rows[i].size = 0;
        continue;
      }
      whichRowAddsToRowsBelow = i;
      dT->reset(maxCore, dim - i - 1);
      parallelFor(0, maxCore, *this);
    }
  }
};


// [[Rcpp::export]]
NumericVector solveSparseMat(NumericMatrix Xmat, NumericVector y, int maxCore = 1)
{
  List tmp = fullMatToSparseColRowContent(Xmat);


  IntegerVector rowSizes = tmp[0];
  NumericVector rowContents = tmp[1];
  INT nrow = y.size();
  vec<valtype> ycopy(y.begin(), y.end());
  vec<sparseV> sparseM(nrow);
  void *rc = &rowContents[0];
  sparseM[0].size = rowSizes[0];
  sparseM[0].ind = adrs<indtype> (rc);
  sparseM[0].val = adrs<valtype> (sparseM[0].ind + sparseM[0].size);
  for(INT i = 1; i < nrow; ++i)
  {
    sparseM[i].size = rowSizes[i];
    sparseM[i].ind = adrs<indtype> (
      sparseM[i - 1].val + sparseM[i - 1].size);
    sparseM[i].val = adrs<valtype> (
      sparseM[i].ind + sparseM[i].size);
  }


  for(int i = 0; i < (int)nrow; ++i)
  {
    for(int j = 0; j < sparseM[i].size; ++j)
    {
      Rcout << sparseM[i].ind[j] << " ";
    }
    Rcout << " ---- ";
    for(int j = 0; j < sparseM[i].size; ++j)
    {
      Rcout << sparseM[i].val[j] << " ";
    }
    Rcout << "\n";
  }


  Rcout << "==================================================================\n";


  valtype *Y = &ycopy[0];
  NumericVector solution(nrow);
  sparseGauEli(nrow, &sparseM[0], Y, &solution[0], maxCore);


  for(int i = 0; i < (int)nrow; ++i)
  {
    for(int j = 0; j < sparseM[i].size; ++j)
    {
      Rcout << sparseM[i].ind[j] << " ";
    }
    Rcout << " ---- ";
    for(int j = 0; j < sparseM[i].size; ++j)
    {
      Rcout << sparseM[i].val[j] << " ";
    }
    Rcout << "\n";
  }


  // return 0;


  for(indtype i = nrow - 1; i >= 0; --i)
  {
    sparseV &row = sparseM[i];
    if(row.size == 0) continue;
    if(row.size == 1)
    {
      solution[row.ind[0]] = Y[i] / row.val[0];
      continue;
    }
    valtype S = 0;
    for(indtype j = 1; j < row.size; ++j)
    {
      S += solution[row.ind[j]] * row.val[j];
    }
    solution[row.ind[0]] = (Y[i] - S) / row.val[0];
  }
  return solution;
}




// ============================================================================
// END:
// Parallel Gaussian elimination to triangle matrix U then back-substitution
// for sparse matrix
// ============================================================================
*/








// Gauss-seidel
valtype innerProd(sparseV &x, valtype *y)
{
  valtype S = 0;
  for(indtype i = 0; i < x.size; ++i)
  {
    S += x.val[i] * y[x.ind[i]];
  }
  return S;
}




// Return true if converged.
// y will change
bool gaussSeidelSparse(valtype *rst, sparseV *rows, valtype *y, indtype dim,
                 valtype convergenceEPS, indtype maxIter,
                 indtype &iter, valtype &maxErr)
{
  // vec<indtype> map(dim);
  vec<unsigned char> mapped(dim, false);
  for(indtype i = 0; i < dim; ++i)
  {
    sparseV &x = rows[i];
    for(indtype j = 0; j < x.size; ++j)
    {
      if(mapped[x.ind[j]]) continue;
      mapped[x.ind[j]] = true;
      // map[i] = x.ind[j];
      valtype r = -1.0 / x.val[j];
      for(indtype k = 0; k < x.size; ++k)
      {
        x.val[k] *= r;
      }
      std::copy(x.ind + j + 1, x.ind + x.size, x.ind + j);
      std::copy(x.val + j + 1, x.val + x.size, x.val + j);
      --x.size;
      y[i] *= r;
      break;
    }
  }


  std::copy(y, y + dim, rst);
  iter = 0;
  while(iter < maxIter)
  {
    maxErr = 0;
    for(indtype k = 0; k < dim; ++k)
    {
      valtype pre = rst[k];
      rst[k] = innerProd(rows[k], rst) - y[k];
      maxErr = std::max<valtype> (maxErr, std::abs(relaErr(pre, rst[k])));
    }
    if(maxErr < convergenceEPS) break;
    ++iter;
  }


  return maxErr < convergenceEPS;
}


// [[Rcpp::export]]
List fullMatToSparseRow(NumericMatrix Xmat)
{
  INT ncol = Xmat.ncol(), nrow = Xmat.nrow();
  // INT nonZero = 0;
  double *X = &Xmat[0];
  IntegerVector sizes(nrow);
  INT containerSize = 0;
  for(INT i = 0; i < nrow; ++i)
  {
    INT rowNonZero = 0;
    for(INT j = 0; j < ncol; ++j)
    {
      if(X[j * nrow + i] == 0) continue;
      ++rowNonZero;
    }
    sizes[i] = rowNonZero;
    INT n = (sizeof(indtype) + sizeof(valtype)) * rowNonZero;
    INT q = n / sizeof(double), r = n % sizeof(double);
    if(r != 0) ++q;
    containerSize += q;
  }
  NumericVector content(containerSize);
  void *rst = &content[0];
  for(INT i = 0; i < nrow; ++i)
  {
    indtype *ind = adrs<indtype> (rst);
    INT k = 0;
    for(INT j = 0; j < ncol; ++j)
    {
      if(X[j * nrow + i] == 0) continue;
      ind[k] = j;
      ++k;
    }
    valtype *val = adrs<valtype> (ind + sizes[i]);
    for(INT j = 0, jend = sizes[i]; j < jend; ++j)
    {
      val[j] = X[ind[j] * nrow + i];
    }
    rst = (INT*)(val + sizes[i]);
  }
  return List::create(Named("sizes") = sizes, Named("content") = content);
}


// [[Rcpp::export]]
List fullMatToSparseRowArmaFormat(arma::mat &Xmat)
{
  INT ncol = Xmat.n_cols, nrow = Xmat.n_rows;
  double *X = &Xmat[0];
  IntegerVector sizes(nrow);
  INT containerSize = 0;
  for(INT i = 0; i < nrow; ++i)
  {
    INT rowNonZero = 0;
    for(INT j = 0; j < ncol; ++j)
    {
      if(X[j * nrow + i] == 0) continue;
      ++rowNonZero;
    }
    sizes[i] = rowNonZero;
    INT n = (sizeof(indtype) + sizeof(valtype)) * rowNonZero;
    INT q = n / sizeof(double), r = n % sizeof(double);
    if(r != 0) ++q;
    containerSize += q;
  }
  NumericVector content(containerSize);
  void *rst = &content[0];
  for(INT i = 0; i < nrow; ++i)
  {
    indtype *ind = adrs<indtype> (rst);
    INT k = 0;
    for(INT j = 0; j < ncol; ++j)
    {
      if(X[j * nrow + i] == 0) continue;
      ind[k] = j;
      ++k;
    }
    valtype *val = adrs<valtype> (ind + sizes[i]);
    for(INT j = 0, jend = sizes[i]; j < jend; ++j)
    {
      val[j] = X[ind[j] * nrow + i];
    }
    rst = (INT*)(val + sizes[i]);
  }
  return List::create(Named("sizes") = sizes, Named("content") = content);
}


// [[Rcpp::export]]
NumericVector guassSeidelSolveSparseDenseInput(
    NumericMatrix X, NumericVector y, int maxIter = 1e6, double convergenceEPS = 1e-6)
{
  List tmp = fullMatToSparseRow(X);
  IntegerVector rowSizes = tmp[0];
  NumericVector rowContents = tmp[1];
  INT nrow = y.size();
  vec<sparseV> sparseM(nrow);
  void *rc = &rowContents[0];
  sparseM[0].size = rowSizes[0];
  sparseM[0].ind = adrs<indtype> (rc);
  sparseM[0].val = adrs<valtype> (sparseM[0].ind + sparseM[0].size);
  for(INT i = 1; i < nrow; ++i)
  {
    sparseM[i].size = rowSizes[i];
    sparseM[i].ind = adrs<indtype> (
      sparseM[i - 1].val + sparseM[i - 1].size);
    sparseM[i].val = adrs<valtype> (
      sparseM[i].ind + sparseM[i].size);
  }


  NumericVector rst(y.size());
  indtype iter = 0;
  valtype maxErr = 0;
  bool converged = gaussSeidelSparse(
    &rst[0], &sparseM[0], &(vec<valtype> (y.begin(), y.end())[0]),
    nrow, convergenceEPS, maxIter, iter, maxErr);
  Rcout << "Iteration = " << iter;
  if(converged) Rcout << ", converged.\n";
  else Rcout << ", unconverged! Maximal error = " << maxErr << "\n";
  return rst;
}


arma::mat guassSeidelSolveSparseDenseInputArmaFormat(
    arma::mat &X, arma::colvec &y, int maxIter = 1e6, double convergenceEPS = 1e-6)
{
  List tmp = fullMatToSparseRowArmaFormat(X);
  IntegerVector rowSizes = tmp[0];
  NumericVector rowContents = tmp[1];
  INT nrow = y.size();
  vec<sparseV> sparseM(nrow);
  void *rc = &rowContents[0];
  sparseM[0].size = rowSizes[0];
  sparseM[0].ind = adrs<indtype> (rc);
  sparseM[0].val = adrs<valtype> (sparseM[0].ind + sparseM[0].size);
  for(INT i = 1; i < nrow; ++i)
  {
    sparseM[i].size = rowSizes[i];
    sparseM[i].ind = adrs<indtype> (
      sparseM[i - 1].val + sparseM[i - 1].size);
    sparseM[i].val = adrs<valtype> (
      sparseM[i].ind + sparseM[i].size);
  }


  arma::mat rst(y.size(), 1);
  indtype iter = 0;
  valtype maxErr = 0;
  bool converged = gaussSeidelSparse(
    &rst[0], &sparseM[0], &(vec<valtype> (y.begin(), y.end())[0]),
    nrow, convergenceEPS, maxIter, iter, maxErr);
  Rcout << "Iteration = " << iter;
  if(converged) Rcout << ", converged.\n";
  else Rcout << ", unconverged! Maximal error = " << maxErr << "\n";
  return rst;
}


// [[Rcpp::export]]
NumericVector guassSeidelSolveSparseSparseInput(
    List sparseX, NumericVector y, int maxIter = 1e6, double convergenceEPS = 1e-6)
{
  IntegerVector rowSizes = sparseX[0];
  NumericVector rowContents = sparseX[1];
  INT nrow = y.size();
  vec<sparseV> sparseM(nrow);
  void *rc = &rowContents[0];
  sparseM[0].size = rowSizes[0];
  sparseM[0].ind = adrs<indtype> (rc);
  sparseM[0].val = adrs<valtype> (sparseM[0].ind + sparseM[0].size);
  for(INT i = 1; i < nrow; ++i)
  {
    sparseM[i].size = rowSizes[i];
    sparseM[i].ind = adrs<indtype> (
      sparseM[i - 1].val + sparseM[i - 1].size);
    sparseM[i].val = adrs<valtype> (
      sparseM[i].ind + sparseM[i].size);
  }


  NumericVector rst(y.size());
  indtype iter = 0;
  valtype maxErr = 0;
  bool converged = gaussSeidelSparse(
    &rst[0], &sparseM[0], &(vec<valtype> (y.begin(), y.end())[0]),
    nrow, convergenceEPS, maxIter, iter, maxErr);
  Rcout << "Iteration = " << iter;
  if(converged) Rcout << ", converged.\n";
  else Rcout << ", unconverged! Maximal error = " << maxErr << ".\n";
  return rst;
}








// [[Rcpp::export]]
List paraOLSsparse(NumericVector y, NumericMatrix X,
                   List sparseMat = R_NilValue,
                   double sparseNonZeroThreshold = 0.2,
                   int maxCore = 7, int verboseInt = 8192,
                   bool useGaussSeidel = true, int maxIter = 1e6,
                   double convergenceEPS = 1e-6)
{
  if(sparseMat.size() == 0)
  {
    sparseMat = fullMatToSparseCol(X);
  }


  IntegerVector siz = sparseMat[0];
  IntegerVector ind = sparseMat[1];
  NumericVector val = sparseMat[2];
  INT ncol = siz.size();
  INT nrow = y.size();
  vec<sparseV> sparseM(ncol);


  sparseM[0].size = siz[0];
  sparseM[0].ind = &ind[0];
  sparseM[0].val = &val[0];
  for(INT i = 1; i < ncol; ++i)
  {
    sparseM[i].size = siz[i];
    sparseM[i].ind = sparseM[i - 1].ind + sparseM[i - 1].size;
    sparseM[i].val = sparseM[i - 1].val + sparseM[i - 1].size;
  }


  arma::mat XtXmat(ncol, ncol);
  XtXsparse(verboseInt, ncol, &sparseM[0], &*XtXmat.begin(), maxCore);


  arma::colvec XtYvec(ncol);
  XtY(nrow, ncol, &*X.begin(), &*y.begin(), &*XtYvec.begin(), maxCore);
  if(verboseInt > 0) Rcout << "Solving symmetric linear system...\n";
  // Check if the inner-product square matrix is sparse
  bool isSparse = false;
  {
    INT Nnonzero = 0;
    valtype *x = &*XtXmat.begin();
    INT size = ncol * ncol;
    for(INT i = 0; i < size; ++i)
    {
      if(x[i] != 0) ++Nnonzero;
    }
    if((Nnonzero + 0.0) / size <= sparseNonZeroThreshold) isSparse = true;
    Rcout << "Nonzero entries constitute " << (Nnonzero + 0.0) / size * 100 << "% of the Gram matrix\n";
  }
  arma::mat coef;
  if(!!useGaussSeidel)
  {
    arma::solve(coef, XtXmat, XtYvec);
  }
  else
  {
    if(!isSparse)
    {
      arma::mat tmp = guassSeidelSolveDense(XtXmat, XtYvec, maxIter, convergenceEPS);
      coef.swap(tmp);
    }
    else
    {
      arma::mat tmp = guassSeidelSolveSparseDenseInputArmaFormat(
        XtXmat, XtYvec, maxIter, convergenceEPS);
      coef.swap(tmp);
    }
  }


  NumericVector pred(nrow);
  calPred(nrow, ncol, &*X.begin(), &*coef.begin(), &*pred.begin(), maxCore);
  double totalErr = 0;
  for(INT i = 0; i < nrow; ++i)
  {
    totalErr += (pred[i] - y[i]) * (pred[i] - y[i]);
  }
  return List::create(Named("coef") = coef, Named("pred") = pred,
                      Named("totalErr") = totalErr);
}




































