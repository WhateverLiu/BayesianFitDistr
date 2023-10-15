// #pragma once


// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;


#ifndef losstype
#define losstype double
#endif


#ifndef probtype
#define probtype double
#endif


#ifndef losstype
#define losstype double
#endif


#ifndef vec
#define vec std::vector
#endif


#ifndef EPS100K
#define EPS100K 1e-12
#endif


// #include "h/checkPMFinput.hpp"
extern void checkPMFintegrity(Rcpp::List X, bool checkReverseOrder = false);
extern void checkPMFlistIntegrity(Rcpp::List distlist, bool checkReverseOrder = false);



bool rg4details(losstype *inLoss, probtype *inProb, int inSize, losstype *outLoss,
         probtype *outProb, int outSize,
         int &countOfNegReDist) 
{
  
  bool unableToWork = 0;
  countOfNegReDist = 0;
  
  std::vector<probtype> safe(outProb, outProb + outSize);
  
  if (outSize < 5 || inSize < 5) 
  {
    unableToWork = 1;
    // DownSample(inLoss, inProb, inSize, outLoss, outProb, outSize);
  } 
  else
  {
    // std::fill(outProb, outProb + outSize, (double)0);
    
    int n = 0;
    
    //// x3 is third loss points. Either last or first loss
    // Loss_t x3 = outLoss[outSize - 1];
    //// up to the last interval
    losstype x0 = outLoss[0];
    losstype x3 = outLoss[outSize - 1];
    
    probtype p0 = 0;
    probtype p1 = 0;
    probtype p2 = 0;
    probtype p3 = 0;
    probtype p = 0;
    
    losstype x = 0;
    losstype x1 = 0;
    losstype x2 = 0;
    for (int i = 0, iend = outSize - 1; i < iend;) {
      if (std::abs(inLoss[n] - outLoss[i]) < EPS100K) {
        outProb[i] += inProb[n];
        ++n;
        if (n >= inSize) {
          break;
        }
      } else {
        x1 = outLoss[i];
        x2 = outLoss[i + 1];
        
        
        while (n < inSize && (x2 - inLoss[n]) > EPS100K)
        {
          x = inLoss[n];
          p = inProb[n];
          double tmp = p * (x3 - x) * (x0 - x) /
            ((x3 - x1) * (x0 - x1) * (x2 - x) +
              (x3 - x2) * (x0 - x2) * (x - x1));
          p1 = tmp * (x2 - x);
          p2 = tmp * (x - x1);
          p0 = (p * (x3 - x) - (x3 - x1) * p1 - (x3 - x2) * p2) / (x3 - x0);
          p3 = p - p0 - p1 - p2;
          ++n;
          outProb[0] += p0;
          outProb[i] += p1;
          outProb[i + 1] += p2;
          outProb[outSize - 1] += p3;
        }
        i++;
      }
    }
    
    
    if (std::abs(outLoss[outSize - 1] - inLoss[inSize - 1]) < EPS100K)
      outProb[outSize - 1] += inProb[inSize - 1];
    
    
    if (outProb[0] >= 0 && outProb[outSize - 1] >= 0) return 0;
    
    
    // bool unableToWork = false;
    
    int i = 0;
    int j = outSize - 1;
    
    
    while (true) 
    {
      
      // for (int u = 0; u < outSize; ++u)
      //   Rcout << outProb[u] << ", ";
      // Rcout << "\n\n";
      
      
      if (outProb[i] < 0) {
        
        countOfNegReDist += 1;
        x = outLoss[i];
        x1 = outLoss[i + 1];
        x2 = outLoss[i + 2];
        x3 = outLoss[j];
        
        p = outProb[i];
        
        p1 = (x - x3) * (x - x2) / ((x1 - x3) * (x1 - x2)) * p;
        p2 = (x - x1) * (x - x3) / ((x2 - x3) * (x2 - x1)) * p;
        p3 = p - p1 - p2;
        outProb[i + 1] += p1;
        outProb[i + 2] += p2;
        outProb[j] += p3;
        outProb[i] = 0;
        ++i;
      }
      
      
      if (outProb[i] >= 0 && outProb[j] >= 0) {
        break;
      }
      
      
      if (j - i < 3) {
        unableToWork = true;
        break;
      }
      
      
      if (outProb[j] < 0) {
        countOfNegReDist += 1;
        
        
        x = outLoss[j];
        x1 = outLoss[j - 1];
        x2 = outLoss[j - 2];
        x3 = outLoss[i];
        p = outProb[j];
        
        
        p1 = (x - x3) * (x - x2) / ((x1 - x3) * (x1 - x2)) * p;
        p2 = (x - x1) * (x - x3) / ((x2 - x3) * (x2 - x1)) * p;
        p3 = p - p1 - p2;
        outProb[j - 1] += p1;
        outProb[j - 2] += p2;
        outProb[i] += p3;
        outProb[j] = 0;
        --j;
      }
      
      
      if (outProb[i] >= 0 && outProb[j] >= 0) {
        break;
      }
      
      
      if (j - i < 3) {
        unableToWork = true;
        break;
      }
    }
  }
  
  
  if (unableToWork) 
  {
    int j = 0;
    
    std::copy(safe.begin(), safe.end(), outProb);
    
    double tmp, interval;
    
    losstype *ng_st = outLoss, *x_st = inLoss;
    
    probtype *ngp_st = outProb, *xp_st = inProb;
    
    int ngsize = outSize, sizex = inSize;
    
    for (int i = 1; i < ngsize; ++i) 
    {
      interval = *(ng_st + i) - *(ng_st + i - 1);
      while (j != sizex && *(x_st + j) <= *(ng_st + i)) {
        tmp = *(xp_st + j) / interval * (*(ng_st + i) - *(x_st + j));
        *(ngp_st + i - 1) += tmp;
        *(ngp_st + i) += *(xp_st + j) - tmp;
        ++j;
      }
    }
  }
  
  
  return unableToWork;
}


//'
//' Four-point regriding details
//' 
//' Regrid a PMF using four-point regriding. The algorithm usually but not always 
//' preserves the PMF's second moment. In addition to the regrided PMF, the function 
//' returns some process information.
//' 
//' @inherit rglr4split
//' 
//' @return A list of 3:
//' \itemize{
//'   \item{\code{$RTL}}{ TRUE if the algorithm retreated to linear regriding.}
//'   \item{\code{$Nredist}}{ An integer. The number of redistributions of negative 
//'   probabilities in Stage II.}
//'   \item{\code{$resultDistr}}{ A 2-column data frame 
//'   \code{(support, probability)} as the regrided PMF.}
//'  }
//' 
//' @example inst/examples/rglr4splitDetails.R
// [[Rcpp::export]]
List rglr4splitDetails(DataFrame X, NumericVector ngd) 
{
  checkPMFintegrity(X);
  NumericVector val = X[0], P = X[1], newgridP(ngd.size());
  int negCount = 0;
  bool retreatedToLinearRegrid = rg4details(
    &val[0], &P[0], val.size(), &ngd[0], 
    &newgridP[0], ngd.size(), negCount);
  DataFrame rstDistr = DataFrame::create(
    Named("val") = ngd, Named("P") = newgridP);
  return List::create(Named("RTL") = retreatedToLinearRegrid,
                      Named("Nredist") = negCount,
                      Named("resultDistr") = rstDistr);
}




#ifdef losstype
#undef losstype
#endif


#ifdef probtype
#undef probtype
#endif


#ifdef vec
#undef vec
#endif


#ifdef EPS100K
#undef EPS100K
#endif














