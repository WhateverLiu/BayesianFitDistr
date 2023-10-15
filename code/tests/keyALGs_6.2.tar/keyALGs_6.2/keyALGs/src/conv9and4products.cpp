// #pragma once

// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>

#include <cmath>
#include <iostream>
// #include <fstream>
#include <RcppParallel.h>
#include <setjmp.h>
#include <string.h>

jmp_buf env;

// [[Rcpp::depends(RcppParallel)]]
using namespace RcppParallel;
using namespace Rcpp;

#define EPS100K 1e-12
//#define double double

#define losstype double
#define probtype double
#define vec std::vector
#define indtype unsigned

// unsigned lrCount=0;

// #include "h/checkPMFinput.hpp"
extern void checkPMFintegrity(Rcpp::List X, bool checkReverseOrder = false);
extern void checkPMFlistIntegrity(Rcpp::List distlist, bool checkReverseOrder = false);
  


bool rglinearAddCppPtr(losstype *x_st,  probtype *xp_st,  unsigned sizex,
                       losstype *ng_st, probtype *ngp_st, unsigned ngsize)
{
  if (*(ng_st + ngsize - 1) < *(x_st + sizex - 1))
    *(ng_st + ngsize - 1) = *(x_st + sizex - 1);
  unsigned j = 0;
  losstype tmp, interval /*=(*(ng_st+ngsize-1)-*ng_st)/(ngsize-1)*/;
  for (unsigned i = 1; i != ngsize; ++i) {
    interval = *(ng_st + i) - *(ng_st + i - 1);
    while (*(x_st + j) <= *(ng_st + i) && j != sizex) {
      tmp = *(xp_st + j) / interval * (*(ng_st + i) - *(x_st + j));
      *(ngp_st + i - 1) += tmp;
      *(ngp_st + i) += *(xp_st + j) - tmp;
      ++j;
    }
  }
  return 1;
}


bool rglinearAddCpp(std::vector<losstype>::iterator x_st,
                    std::vector<probtype>::iterator xp_st, unsigned sizex,
                    std::vector<losstype>::iterator ng_st,
                    std::vector<probtype>::iterator ngp_st, unsigned ngsize)
// warning!: this function will check ng_st.back() and *(x_st+sizex-1)
// immediately if ng_st.back() is less, it will set it to *(x_st+sizex-1)!
{
  return rglinearAddCppPtr(&*x_st, &*xp_st, sizex, &*ng_st, &*ngp_st, ngsize);
}


bool rglr3splitAddCpp(std::vector<double>::iterator x_st,
                      std::vector<double>::iterator xp_st, unsigned sizex,
                      std::vector<double>::iterator ng_st,
                      std::vector<double>::iterator ngp_st_final,
                      unsigned ngsize)
// warning!: this function will check ng_st.back() and *(x_st+sizex-1)
// immediately if ng_st.back() is less, it will set it to *(x_st+sizex-1)!
{
  if (*(ng_st + ngsize - 1) < *(x_st + sizex - 1))
    *(ng_st + ngsize - 1) = *(x_st + sizex - 1);
  std::vector<double> ngpTmp(ngsize, 0);
  std::vector<double>::iterator ngp_st = ngpTmp.begin();
  unsigned j = 0;
  std::vector<double>::iterator first = ng_st, last = ng_st + ngsize - 1;
  double middle = (*first + *last) / 2;
  unsigned firstHalfEndX = std::lower_bound(x_st, x_st + sizex, middle) - x_st;
  std::vector<double>::iterator endOfNgp = ngp_st + ngsize - 1,
                                endOfNg = ng_st + ngsize - 1, farGrid;
  std::vector<double>::iterator ng_stReserve = ng_st /*, ngp_stReserve=ngp_st*/;
  double x1, x2, x3, x, p, p1, p2, p3;
  bool indicator = 0;
  if (ngsize < 4 || sizex < 4) indicator = 1;
  if (!indicator) {
    for (unsigned i = 1; i != ngsize; ++i) {
      x1 = *(ng_st + i - 1);
      x2 = *(ng_st + i);
      while (*(x_st + j) <= *(ng_st + i) && j != sizex) {
        if (j < firstHalfEndX) {
          x3 = *endOfNg;
          farGrid = endOfNgp;
        } else {
          x3 = *ng_st;
          farGrid = ngp_st;
        }
        x = *(x_st + j);
        p = *(xp_st + j);
        p1 = (x - x3) * (x - x2) / ((x1 - x3) * (x1 - x2)) * p;
        p2 = (x - x1) * (x - x3) / ((x2 - x3) * (x2 - x1)) * p;
        p3 = p - p1 - p2;
        *(ngp_st + i - 1) += p1;
        *(ngp_st + i) += p2;
        *farGrid += p3;
        ++j;
      }
    }
    while (1) {
      if (*ngp_st < 0) {
        x = *ng_st;
        x1 = *(ng_st + 1);
        x2 = *(ng_st + 2);
        x3 = *endOfNg;
        p = *ngp_st;
        p1 = (x - x3) * (x - x2) / ((x1 - x3) * (x1 - x2)) * p;
        p2 = (x - x1) * (x - x3) / ((x2 - x3) * (x2 - x1)) * p;
        p3 = p - p1 - p2;
        *(ngp_st + 1) += p1;
        *(ngp_st + 2) += p2;
        *endOfNgp += p3;
        ++ng_st;
        *ngp_st = 0;
        ++ngp_st;
      }
      if (*ngp_st >= 0 && *endOfNgp >= 0) break;
      if (endOfNg - ng_st < 3) {
        indicator = 1;
        break;
      }
      if (*endOfNgp < 0) {
        x = *endOfNg;
        x1 = *(endOfNg - 1);
        x2 = *(endOfNg - 2);
        x3 = *ng_st;
        p = *endOfNgp;
        p1 = (x - x3) * (x - x2) / ((x1 - x3) * (x1 - x2)) * p;
        p2 = (x - x1) * (x - x3) / ((x2 - x3) * (x2 - x1)) * p;
        p3 = p - p1 - p2;
        *(endOfNgp - 1) += p1;
        *(endOfNgp - 2) += p2;
        *ngp_st += p3;
        --endOfNg;
        *endOfNgp = 0;
        --endOfNgp;
      }
      if (*ngp_st >= 0 && *endOfNgp >= 0) break;
      if (endOfNg - ng_st < 3) {
        indicator = 1;
        break;
      }
    }
  }
  if (indicator) {
    unsigned j = 0;
    ng_st = ng_stReserve;
    ngp_st = ngp_st_final;
    double tmp, interval;
    for (unsigned i = 1; i != ngsize; ++i) {
      interval = *(ng_st + i) - *(ng_st + i - 1);
      while (*(x_st + j) <= *(ng_st + i) && j != sizex) {
        tmp = *(xp_st + j) / interval * (*(ng_st + i) - *(x_st + j));
        *(ngp_st + i - 1) += tmp;
        *(ngp_st + i) += *(xp_st + j) - tmp;
        ++j;
      }
    }
  } else {
    ++endOfNgp;
    for (ngp_st_final += ngp_st - ngpTmp.begin(); ngp_st != endOfNgp;
         ++ngp_st, ++ngp_st_final)
      *ngp_st_final += *ngp_st;
  }
  return indicator;
}




// x_st: the old distribution's variable container head pointer.
// xp_st: the old distribution's probability container head pointer.
// sizex: length of the old distribution.
// ng_st: the new grid's container head pointer. This container should have
// contained the new distribution's variable values ngp_st: the new probability
// container head pointer. This container should have container 0s. ngsize:
// length of the new distribution.
bool rg4(losstype *inLoss, probtype *inProb, int inSize, losstype *outLoss,
         probtype *outProb, int outSize) 
{

  bool unableToWork = 0;

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
    for (int i = 0; i < outSize - 1;) {
      if (std::abs(inLoss[n] - outLoss[i]) < EPS100K) {
        // if (inLoss[n]==outLoss[i]) {
        outProb[i] += inProb[n];
        ++n;
        if (n >= inSize) {
          break;
        }
      } else {
        x1 = outLoss[i];
        x2 = outLoss[i + 1];

        ////while (n < inSize && inLoss[n] < x2 - EPS100K)
        while (n < inSize && (x2 - inLoss[n]) > EPS100K)
        // while (n < inSize && x2>=inLoss[n])
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
      if (outProb[i] < 0) {
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

    for (int i = 1; i != ngsize; ++i) 
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




bool rglr4splitAddCpp(std::vector<losstype>::iterator x_st,
                      std::vector<probtype>::iterator xp_st, unsigned sizex,
                      std::vector<losstype>::iterator ng_st,
                      std::vector<probtype>::iterator ngp_st, unsigned ngsize) 
{
  return rg4(&*x_st, &*xp_st, sizex, &*ng_st, &*ngp_st, ngsize);
}




bool lmmCpp(losstype *x, probtype *p, int siz, losstype *xnew, probtype *pnew,
            int sizeNew) 
{
  if (xnew[sizeNew - 1] < x[siz - 1]) xnew[sizeNew - 1] = x[siz - 1];
  if (xnew[0] > x[0]) xnew[0] = x[0];

  losstype *left = &xnew[0], *xnewEnd = xnew + sizeNew;
  probtype *leftP = &pnew[0];
  // mid=left[1], right=left[2], rightright=left[3];

  int i = 0, iend = siz;  // i points to x[i]

  losstype x1_x0 = left[1] - left[0], x2_x1 = left[2] - left[1],
           x2_x0 = x2_x1 + x1_x0, leftHeadTailMid = (left[0] + left[3]) / 2;

  vec<probtype> safe(pnew, pnew + sizeNew);

  bool notWorking = 0;

  if (siz < 4 or sizeNew < 4)
    notWorking = 1;
  else {
    for (; i < iend; ++i) {
      if (x[i] > leftHeadTailMid)  // x[i]-left[0]>left[3]-x[i]
      {
        if (leftP[0] < 0) {
          notWorking = 1;
          break;
        }
        ++left;

        if (left + 2 >= xnewEnd) {
          --left;
          break;
        }

        ++leftP;

        leftHeadTailMid = (left[0] + left[3]) / 2;
        x1_x0 = x2_x1;
        x2_x1 = left[2] - left[1];
        x2_x0 = x2_x1 + x1_x0;
      }

      leftP[0] += (x[i] - left[1]) * (x[i] - left[2]) / (x1_x0 * x2_x0) * p[i];
      leftP[1] += (x[i] - left[0]) * (x[i] - left[2]) / (-x1_x0 * x2_x1) * p[i];
      leftP[2] += (x[i] - left[0]) * (x[i] - left[1]) / (x2_x0 * x2_x1) * p[i];
    }
  }

  if (!notWorking and i != iend) {
    for (; i < iend; ++i) {
      leftP[0] += (x[i] - left[1]) * (x[i] - left[2]) / (x1_x0 * x2_x0) * p[i];
      leftP[1] += (x[i] - left[0]) * (x[i] - left[2]) / (-x1_x0 * x2_x1) * p[i];
      leftP[2] += (x[i] - left[0]) * (x[i] - left[1]) / (x2_x0 * x2_x1) * p[i];
    }
  }

  if (!notWorking && leftP[0] >= 0 && leftP[1] >= 0 && leftP[2] >= 0) return 1;

  // probability container polluted, copy the safe and start linear regridding

  std::copy(safe.begin(), safe.end(), pnew);


  int j = 0;
  for (int i = 1; i != sizeNew; ++i) 
  {
    losstype interval = xnew[i] - xnew[i - 1];
    while (j != siz && x[j] <= xnew[i]) 
    {
      probtype tmp = p[j] / interval * (xnew[i] - x[j]);
      pnew[i - 1] += tmp;
      pnew[i] += p[j] - tmp;
      ++j;
    }
  }

  return 0;
}




bool lmmAddCpp(std::vector<losstype>::iterator x_st,
               std::vector<probtype>::iterator xp_st, unsigned sizex,
               std::vector<losstype>::iterator ng_st,
               std::vector<probtype>::iterator ngp_st, unsigned ngsize) 
{
  return lmmCpp(&*x_st, &*xp_st, sizex, &*ng_st, &*ngp_st, ngsize);
}




// [[Rcpp::export]]
List lmm2ndM(DataFrame dist, NumericVector newSupport) 
{
  checkPMFintegrity(dist);
  NumericVector val = dist[0], P = dist[1];
  NumericVector newP(newSupport.size());
  losstype *x = &val[0], *xnew = &newSupport[0];
  probtype *p = &P[0], *pnew = &newP[0];
  bool rst = lmmCpp(x, p, val.size(), xnew, pnew, newP.size());
  if (!rst) Rcout << "LMM failed. Please use linear regrid.\n";
  return DataFrame::create(Named("val") = newSupport, Named("P") = newP);
}




// DONT DELETE!!!
// bool rglr4splitAddCpp(std::vector<double>::iterator x_st,
// std::vector<double>::iterator xp_st,
//               unsigned sizex, std::vector<double>::iterator ng_st,
//               std::vector<double>::iterator ngp_st, unsigned ngsize){
//
// // std::ofstream myfile;
// // myfile.open("C:/Users/i56087/Desktop/keyALGs/packages/New
// folder/error.csv");
// // myfile<<"openned!"<<std::endl;
//
//
// // check if the new grid covers the old grid. If not, force it.
// if(*(ng_st+ngsize-1)<*(x_st+sizex-1))*(ng_st+ngsize-1)=*(x_st+sizex-1);
// if(*ng_st>*x_st)*ng_st=*x_st;
//
// std::vector<double>safe(ngp_st, ngp_st+ngsize);
//
// //double safe[ngsize], *ngp_stPointer=&*ngp_st;
//
// //memcpy(safe, ngp_stPointer, sizeof(double)*ngsize);
//
// unsigned i=1, j=0;
// std::vector<double>::iterator endOfNgp=ngp_st+ngsize-1,
// endOfNg=ng_st+ngsize-1; std::vector<double>::iterator ng_stReserve=ng_st,
// ngp_stReserve=ngp_st; double x0=*ng_st, x1, x2, x3=*endOfNg, x, p, p0, p1,
// p2, p3; bool unableToWork=0;
//
// // if distribution's size < 5 or the new grid has < 5 points, new regridding
// alg is guaranteed to fail. if(ngsize<5||sizex<5)unableToWork=1;
//
// if(!unableToWork)
// {
//   unsigned jend=sizex, iend=ngsize;
//
  // The following 2 IF statements are to avoid 0 denominator in the formula
  // for regridding
  // if( std::abs( *ng_st - *x_st )<1e-20 ) // There must be your way to set the minimal error
  // {
  //   *ngp_st += *xp_st;
  //   ++j;
  //   if ( *(x_st + 1) > *(ng_st + 1) ) ++i;
  // }

//   if(std::abs(*(ng_st+ngsize-1)-*(x_st+sizex-1))<1e-20)//there must be your
//   way to set the minimal error
//   {
//     *(ngp_st+ngsize-1)+=*(xp_st+sizex-1);
//     --jend;
//     //if(*(x_st+sizex-2)<=*(ng_st+ngsize-2))--iend;
//   }
//
//   //core regridding procedure
//   for(;i!=iend;++i)
//   {
//     x1=*(ng_st+i-1);
//     x2=*(ng_st+i);
//     while(j!=jend&&*(x_st+j)<=x2)
//     {
//       x=*(x_st+j);
//       p=*(xp_st+j);
//       double
//       tmp=p*(x3-x)*(x0-x)/((x3-x1)*(x0-x1)*(x2-x)+(x3-x2)*(x0-x2)*(x-x1));
//       p1=tmp*(x2-x);
//       p2=tmp*(x-x1);
//       p0=(p*(x3-x)-(x3-x1)*p1-(x3-x2)*p2)/(x3-x0);
//       p3=p-p0-p1-p2;
//       *(ngp_st+i-1)+=p1;
//       *(ngp_st+i)+=p2;
//       *endOfNgp+=p3;
//       *ngp_st+=p0;
//       ++j;
//     }
//   }
//
//   //start negating the negative boundary probabilities if there's any
//   while(true)
//   {
//     if(*ngp_st<0)
//     {
//       x=*ng_st;
//       x1=*(ng_st+1);
//       x2=*(ng_st+2);
//       x3=*endOfNg;
//       p=*ngp_st;
//       p1=(x-x3)*(x-x2)/((x1-x3)*(x1-x2))*p;
//       p2=(x-x1)*(x-x3)/((x2-x3)*(x2-x1))*p;
//       p3=p-p1-p2;
//       *(ngp_st+1)+=p1;
//       *(ngp_st+2)+=p2;
//       *endOfNgp+=p3;
//       ++ng_st;
//       *ngp_st=0;
//       ++ngp_st;
//     }
//     if(*ngp_st>=0&&*endOfNgp>=0)break;
//     if(endOfNg-ng_st<3)
//     {
//       unableToWork=1;
//       break;
//     }
//     if(*endOfNgp<0)
//     {
//       x=*endOfNg;
//       x1=*(endOfNg-1);
//       x2=*(endOfNg-2);
//       x3=*ng_st;
//       p=*endOfNgp;
//       p1=(x-x3)*(x-x2)/((x1-x3)*(x1-x2))*p;
//       p2=(x-x1)*(x-x3)/((x2-x3)*(x2-x1))*p;
//       p3=p-p1-p2;
//       *(endOfNgp-1)+=p1;
//       *(endOfNgp-2)+=p2;
//       *ngp_st+=p3;
//       --endOfNg;
//       *endOfNgp=0;
//       --endOfNgp;
//     }
//     if(*ngp_st>=0&&*endOfNgp>=0)break;
//     if(endOfNg-ng_st<3)
//     {
//       unableToWork=1;
//       break;
//     }
//   }
// }
//
// //if the negative probibilites can't be balancing out, run linear regridding
// if(unableToWork)
// {
//   unsigned j=0;
//   ng_st=ng_stReserve;
//   ngp_st=ngp_stReserve;
//
//   //fill the probability container with 0
//   //for(std::vector<double>::iterator
//   k=ngp_st,end=ngp_st+ngsize;k!=end;++k)*k=0.0;
//
//   std::copy(safe.begin(), safe.end(), ngp_st);
//
//   //memcpy(ngp_stPointer, safe, sizeof(double)*ngsize);
//
//   double tmp, interval;
//   for(unsigned i=1;i!=ngsize;++i)
//   {
//     interval=*(ng_st+i)-*(ng_st+i-1);
//     while(*(x_st+j)<=*(ng_st+i)&&j!=sizex)
//     {
//       tmp=*(xp_st+j)/interval*(*(ng_st+i)-*(x_st+j));
//       *(ngp_st+i-1)+=tmp;
//       *(ngp_st+i)+=*(xp_st+j)-tmp;
//       ++j;
//     }
//   }
// }
// return unableToWork;
// }




bool rglr4splitGaussianKernalAddCpp(std::vector<double>::iterator x_st,
                                    std::vector<double>::iterator xp_st,
                                    unsigned sizex,
                                    std::vector<double>::iterator ng_st,
                                    std::vector<double>::iterator ngp_st_final,
                                    unsigned ngsize)
// warning!: this function will check ng_st.back() and *(x_st+sizex-1)
// immediately if ng_st.back() is less, it will set it to *(x_st+sizex-1)!
{
  // std::cout<<"enter in the regridding function\n";
  if (*(ng_st + ngsize - 1) < *(x_st + sizex - 1))
    *(ng_st + ngsize - 1) = *(x_st + sizex - 1);
  std::vector<double> ngpTmp(ngsize, 0);
  std::vector<double>::iterator ngp_st = ngpTmp.begin();

  unsigned j = 0;
  //  std::vector<double>::iterator first=ng_st, last=ng_st+ngsize-1;
  //   double middle=(*first+*last)/2;
  //   unsigned firstHalfEndX=std::lower_bound(x_st, x_st+sizex, middle)-x_st;
  std::vector<double>::iterator endOfNgp = ngp_st + ngsize - 1,
                                endOfNg = ng_st + ngsize -
                                          1 /*, farGrid, secondFarGrid*/;
  std::vector<double>::iterator ng_stReserve = ng_st;
  double x0, x1, x2, x3, x, p, p0, p1, p2, p3;
  bool indicator = 0;
  if (ngsize < 5 || sizex < 5) indicator = 1;

  if (!indicator) {
    for (unsigned i = 1; i != ngsize; ++i) {
      x1 = *(ng_st + i - 1);
      x2 = *(ng_st + i);
      while (*(x_st + j) <= *(ng_st + i) && j != sizex) {
        x3 = *endOfNg;
        x0 = *ng_st;
        x = *(x_st + j);
        p = *(xp_st + j);
        if (std::abs(x2 - x) < 1e-20) {
          *(ngp_st + i) += p;
        } else if (std::abs(x1 - x) < 1e-20) {
          *(ngp_st + i - 1) += p;
        } else {
          //         p1=p*(x3-x)*(x0-x)*(x2-x)/((x3-x1)*(x0-x1)*(x2-x)+(x3-x2)*(x0-x2)*(x-x1));
          //         p2=(x-x1)*p1/(x2-x);
          //         p0=(p*(x3-x)-(x3-x1)*p1-(x3-x2)*p2)/(x3-x0);
          //         p3=p-p0-p1-p2;
          //         *(ngp_st+i-1)+=p1;
          //         *(ngp_st+i)+=p2;
          //         *endOfNgp+=p3;
          //         *ngp_st+=p0;

          double sd = (x2 - x1) / 3, var2 = 2 * sd * sd, c = std::exp(-4.5);
          // std::cout<<"linear k="<<(x-x1)/(x2-x)<<std::endl;
          double k = (std::exp(-(x2 - x) * (x2 - x) / var2) - c) /
                     (std::exp(-(x1 - x) * (x1 - x) / var2) - c);
          //         std::cout<<"k="<<k<<std::endl;
          //         std::cout<<x1<<" "<<x<<" "<<x2<<std::endl;
          //         std::cout<<"exp(-(x1-x)*(x1-x)/var2)="<<std::exp(-(x1-x)*(x1-x)/var2)<<std::endl;
          //         std::cout<<"sd="<<sd<<std::endl;
          //         std::cout<<"var2="<<var2<<std::endl;
          //         std::cout<<"-(3*sd-x)*(3*sd-x)/var2="<<-(3*sd-x)*(3*sd-x)/var2<<std::endl;
          //         std::cout<<"c="<<c<<std::endl;
          // double k=std::exp(4.5*(x1+x2-2*x)/(x1-x2));
          // k=(x-x1)/(x2-x);
          p1 = (x - x3) * (x - x0) * p /
               ((x1 - x3) * (x1 - x0) + k * (x2 - x3) * (x2 - x0));
          p2 = k * p1;
          p0 = (p * (x3 - x) - (x3 - x1) * p1 - (x3 - x2) * p2) / (x3 - x0);
          p3 = p - p0 - p1 - p2;
          *(ngp_st + i - 1) += p1;
          *(ngp_st + i) += p2;
          *endOfNgp += p3;
          *ngp_st += p0;
        }
        ++j;
      }
    }
    while (1) {
      if (*ngp_st < 0) {
        x = *ng_st;
        x1 = *(ng_st + 1);
        x2 = *(ng_st + 2);
        x3 = *endOfNg;
        p = *ngp_st;
        p1 = (x - x3) * (x - x2) / ((x1 - x3) * (x1 - x2)) * p;
        p2 = (x - x1) * (x - x3) / ((x2 - x3) * (x2 - x1)) * p;
        p3 = p - p1 - p2;
        *(ngp_st + 1) += p1;
        *(ngp_st + 2) += p2;
        *endOfNgp += p3;
        ++ng_st;
        *ngp_st = 0;
        ++ngp_st;
      }
      if (*ngp_st >= 0 && *endOfNgp >= 0) break;
      if (endOfNg - ng_st < 3) {
        indicator = 1;
        break;
      }
      if (*endOfNgp < 0) {
        x = *endOfNg;
        x1 = *(endOfNg - 1);
        x2 = *(endOfNg - 2);
        x3 = *ng_st;
        p = *endOfNgp;
        p1 = (x - x3) * (x - x2) / ((x1 - x3) * (x1 - x2)) * p;
        p2 = (x - x1) * (x - x3) / ((x2 - x3) * (x2 - x1)) * p;
        p3 = p - p1 - p2;
        *(endOfNgp - 1) += p1;
        *(endOfNgp - 2) += p2;
        *ngp_st += p3;
        --endOfNg;
        *endOfNgp = 0;
        --endOfNgp;
      }
      if (*ngp_st >= 0 && *endOfNgp >= 0) break;
      if (endOfNg - ng_st < 3) {
        indicator = 1;
        break;
      }
    }
  }
  if (indicator) {
    unsigned j = 0;
    ng_st = ng_stReserve;
    ngp_st = ngp_st_final;
    double tmp, interval;
    for (unsigned i = 1; i != ngsize; ++i) {
      interval = *(ng_st + i) - *(ng_st + i - 1);
      while (*(x_st + j) <= *(ng_st + i) && j != sizex) {
        tmp = *(xp_st + j) / interval * (*(ng_st + i) - *(x_st + j));
        *(ngp_st + i - 1) += tmp;
        *(ngp_st + i) += *(xp_st + j) - tmp;
        ++j;
      }
    }
  } else {
    ++endOfNgp;
    for (ngp_st_final += ngp_st - ngpTmp.begin(); ngp_st != endOfNgp;
         ++ngp_st, ++ngp_st_final)
      *ngp_st_final += *ngp_st;
  }
  return indicator;
}




// [[Rcpp::export]]
DataFrame rglr3split(DataFrame X, NumericVector ngd) 
{
  std::vector<double> val = X[0], P = X[1];
  std::vector<double> newgrid(ngd.begin(), ngd.end()), newgridP(ngd.size(), 0);
  rglr3splitAddCpp(val.begin(), P.begin(), val.size(), newgrid.begin(),
                   newgridP.begin(), newgrid.size());
  return DataFrame::create(Named("val", newgrid), Named("P", newgridP));
}





//'
//' Four-point regriding
//' 
//' Regrid a PMF using four-point regriding. The algorithm usually but not always 
//' preserves the PMF's second moment.
//' 
//' @inherit rglr
//' 
//' @example inst/examples/rglr4split.R
// [[Rcpp::export]]
DataFrame rglr4split(DataFrame X, NumericVector ngd) 
{
  checkPMFintegrity(X);
  NumericVector val = X[0], P = X[1], newgridP(ngd.size());
  rg4(&val[0], &P[0], val.size(), &ngd[0], &newgridP[0], ngd.size()); 
  return DataFrame::create(Named("val") = ngd, Named("P") = newgridP);
}












// [[Rcpp::export]]
DataFrame rglr4splitGuassianKernal(DataFrame X, NumericVector ngd) 
{
  std::vector<double> val = X[0], P = X[1];
  std::vector<double> newgrid(ngd.begin(), ngd.end()), newgridP(ngd.size(), 0);
  rglr4splitGaussianKernalAddCpp(val.begin(), P.begin(), val.size(),
                                 newgrid.begin(), newgridP.begin(),
                                 newgrid.size());
  return DataFrame::create(Named("val", newgrid), Named("P", newgridP));
}





typedef bool (*rgMethod)(std::vector<losstype>::iterator,
                         std::vector<probtype>::iterator, unsigned,
                         std::vector<losstype>::iterator,
                         std::vector<probtype>::iterator, unsigned);




inline unsigned lowerInt(double x) 
{
  int y = x;
  if (std::abs(x - y) < 1e-10) --y;
  if (y < 0) return 0;
  return y;
}


void formGrid(std::vector<losstype>::iterator st, unsigned len,
              losstype initialVal, losstype delta) {
  std::vector<losstype>::iterator end = st + len;
  *st = initialVal;
  for (++st; st != end; ++st) *st = *(st - 1) + delta;
}


void formGrid(std::vector<losstype>::iterator st, double initialVal,
              double endVal, double delta, unsigned &len) {
  len = (endVal - initialVal) / delta;
  std::vector<losstype>::iterator end = st + len;
  ++len;
  *st = initialVal;
  for (; st != end; ++st) *st = *(st - 1) + delta;
  *st = endVal;
}


void constMulV(std::vector<probtype>::iterator st, unsigned len,
               std::vector<probtype>::iterator rst, double constant) {
  std::vector<probtype>::iterator end = st + len;
  for (; st != end; ++st, ++rst) *rst = *st * constant;
}


void constAddV(std::vector<losstype>::iterator st, unsigned len,
               std::vector<losstype>::iterator rst, double constant) {
  std::vector<losstype>::iterator end = st + len;
  for (; st != end; ++st, ++rst) *rst = *st + constant;
}


bool truncateDist(std::vector<losstype> &x, std::vector<probtype> &xp,
                  double headTrunc, double tailTrunc) {
  std::vector<probtype>::iterator i = xp.end() - 1;
  bool activated = 0;
  if (*i < tailTrunc) {
    activated = 1;
    for (--i; i >= xp.begin(); --i) {
      *i += *(i + 1);
      if (*i >= tailTrunc) break;
    }
  }
  ++i;  // now i is the end() of that new vector
  std::vector<probtype>::iterator j = xp.begin();
  if (*j < headTrunc) {
    activated = 1;
    for (++j; j != i; ++j) {
      *j += *(j - 1);
      if (*j >= headTrunc) break;
    }
    std::copy(j, i, xp.begin());
    std::copy(x.begin() + (j - xp.begin()), x.begin() + (i - xp.begin()),
              x.begin());
  }
  xp.resize(i - j);
  x.resize(i - j);
  return activated;
}


bool truncateDistSingle(std::vector<losstype> &x, std::vector<probtype> &xp,
                        double headTrunc, double tailTrunc) {
  std::vector<probtype>::iterator i = xp.end() - 1;
  double S = *i;
  bool activated = 0;
  if (*i < tailTrunc) {
    activated = 1;
    for (--i; i >= xp.begin(); --i) {
      S += *i;
      if (*i >= tailTrunc) break;
    }
  }
  *i = S;
  ++i;  // now i is the end() of that new vector
  std::vector<probtype>::iterator j = xp.begin();
  S = *j;
  if (*j < headTrunc) {
    activated = 1;
    for (++j; j != i; ++j) {
      S += *j;
      if (*j >= headTrunc) {
        *j = S;
        break;
      }
    }
    std::copy(j, i, xp.begin());
    std::copy(x.begin() + (j - xp.begin()), x.begin() + (i - xp.begin()),
              x.begin());
  }
  xp.resize(i - j);
  x.resize(i - j);
  return activated;
}


// [[Rcpp::export]]
DataFrame truncateDistSingleTest(DataFrame X, double headTrunc = 0,
                                 double tailTrunc = 0) {
  std::vector<losstype> x = X[0];
  std::vector<probtype> xp = X[1];
  truncateDistSingle(x, xp, headTrunc, tailTrunc);
  return DataFrame::create(Named("val") = x, Named("P") = xp);
}


// [[Rcpp::export]]
DataFrame truncateDistTest(DataFrame X, double headTrunc = 0,
                           double tailTrunc = 0) {
  std::vector<losstype> x = X[0];
  std::vector<probtype> xp = X[1];
  truncateDist(x, xp, headTrunc, tailTrunc);
  return DataFrame::create(Named("val") = x, Named("P") = xp);
}




// FFT work space
struct complex {
  probtype R, I;
  // complex(){real=0;imagine=0;}
  void copy(complex &x) {
    R = x.R;
    I = x.I;
  }
  void copy(complex *x) {
    R = x->R;
    I = x->I;
  }
};


template <typename T>
void index(T len, T *&parent, T *&child) {
  T parentBlocks = 1;
  while (true) {
    if (parentBlocks == len) break;
    T parentBlockSize = len / parentBlocks;
    for (T *parentBlockSt = parent, *parentBlockStEnd = parent + len;
         parentBlockSt != parentBlockStEnd; parentBlockSt += parentBlockSize) {
      T *childBlockSt = child + (parentBlockSt - parent);
      for (T i = 0, iend = parentBlockSize / 2; i != iend; ++i) {
        childBlockSt[i] = parentBlockSt[2 * i];
        childBlockSt[i + iend] = parentBlockSt[2 * i + 1];
      }
    }
    std::swap(parent, child);
    parentBlocks *= 2;
  }
}

void cmul(complex &rst, complex &x, complex &y) {
  rst.R = x.R * y.R - x.I * y.I;
  rst.I = x.I * y.R + x.R * y.I;
}

// len is twice the size of wFFTuse
void generateMultiplier(unsigned len, complex *wToUseFFT, complex *wToUseIFFT) 
{
  unsigned halfSize = len / 2;
  complex root;
  root.R = std::cos(2 * M_PI / len);
  root.I = std::sqrt(1 - root.R * root.R);

  unsigned twiceM = len * 2;
  std::vector<unsigned> tmp(twiceM);
  unsigned *ind = &tmp.front(), *ind2 = ind + len;

  for (unsigned i = 0; i != len; ++i) ind[i] = i;
  index<unsigned>(len, ind, ind2);

  std::vector<complex> allRootsVec(len);
  complex *allRoots = &allRootsVec.front();

  allRoots->R = 1;
  allRoots->I = 0;
  for (unsigned i = 1; i != len; ++i) cmul(allRoots[i], allRoots[i - 1], root);

  wToUseFFT[0].copy(allRoots[0]);
  wToUseIFFT[0].copy(allRoots[0]);
  for (unsigned i = 1; i < halfSize; ++i) {
    wToUseFFT[i].copy(allRoots[ind[2 * i]]);
    wToUseIFFT[i].copy(allRoots[len - ind[2 * i]]);
  }
}


void fftRoutine(complex &rstUp, complex &rstDown, complex &x, complex &w,
                complex &y) {
  double tmpR = w.R * y.R - w.I * y.I, tmpI = w.I * y.R + w.R * y.I;
  rstUp.R = x.R + tmpR;
  rstUp.I = x.I + tmpI;
  rstDown.R = x.R - tmpR;
  rstDown.I = x.I - tmpI;
}


void fftComplex(complex *&val, unsigned coefSize, complex *wToUse,
                complex *&val2) 
{
  unsigned gap = coefSize / 2;
  while (true) {
    if (gap == 0) break;
    unsigned blockSize = gap * 2;
    for (unsigned i = 0, whichBlock = 0; i != coefSize;
         i += blockSize, ++whichBlock) 
    {
      for (unsigned j = i, jend = i + gap; j != jend; ++j) 
      {
        fftRoutine(val2[j], val2[j + gap], val[j], 
                   wToUse[whichBlock], val[j + gap]);
      }
    }
    std::swap(val, val2);
    gap /= 2;
  }
}


void ifftRoutine(complex &rstUp, complex &rstDown, complex &Ax, complex &A_x,
                 complex &multiplier) 
{
  rstUp.R = (Ax.R + A_x.R) / 2;
  rstUp.I = (Ax.I + A_x.I) / 2;
  complex tmp;
  tmp.R = rstUp.R - A_x.R;
  tmp.I = rstUp.I - A_x.I;
  rstDown.R = tmp.R * multiplier.R - tmp.I * multiplier.I;
  rstDown.I = tmp.I * multiplier.R + tmp.R * multiplier.I;
}


void ifftComplex(complex *&val, unsigned coefSize, complex *wToUse,
                 complex *&val2) 
{
  unsigned gap = 1;
  while (true) 
  {
    if (gap == coefSize) break;
    unsigned blockSize = gap * 2;
    for (unsigned i = 0, whichBlock = 0; i != coefSize;
         i += blockSize, ++whichBlock) 
    {
      for (unsigned j = i, jend = i + gap; j != jend; ++j) 
      {
        ifftRoutine(val2[j], val2[j + gap], val[j], 
                    val[j + gap], wToUse[whichBlock]);
      }
    }
    std::swap(val, val2);
    gap *= 2;
  }
}


void fftConv(complex *&x, complex *&y, complex *&buffer, unsigned len,
             complex *wToUseFFT, complex *wToUseIFFT) {
  fftComplex(x, len, wToUseFFT, buffer);
  fftComplex(y, len, wToUseFFT, buffer);
  for (unsigned i = 0; i != len; ++i) cmul(buffer[i], x[i], y[i]);
  ifftComplex(buffer, len, wToUseIFFT, x);
}


struct wToUseFFT 
{
  // unsigned contentSiz;
  std::vector<std::vector<complex> > content;
  std::vector<unsigned> sizes;
  wToUseFFT() {}
  wToUseFFT(unsigned n) {
    unsigned N = 2;
    unsigned contentSiz = 1;
    while (N < n) {
      N *= 2;
      ++contentSiz;
    }
    content.resize(contentSiz);
    sizes.resize(contentSiz);
    sizes.front() = 2;
    content.front().resize(sizes.front());
    for (unsigned i = 1; i != contentSiz; ++i) {
      sizes[i] = 2 * sizes[i - 1];
      content[i].resize(sizes[i]);
    }
    for (unsigned i = 0; i != contentSiz; ++i) {
      generateMultiplier(content[i].size(), &content[i].front(),
                         &content[i].front() + content[i].size() / 2);
    }
  }
  void giveWtoUse(complex *&wToUseFFT, complex *&wToUseIFFT, unsigned n) {
    std::vector<complex> &tmp =
        content[std::lower_bound(sizes.begin(), sizes.end(), n) -
                sizes.begin()];
    wToUseFFT = &tmp.front();
    wToUseIFFT = wToUseFFT + tmp.size() / 2;
  }
};


void fftConvFromReal(probtype *rst, probtype *p1, probtype *p2, unsigned p1len,
                     unsigned p2len, wToUseFFT &w) {
  complex *wToUseFFT, *wToUseIFFT;
  w.giveWtoUse(wToUseFFT, wToUseIFFT, p1len + p2len - 1);
  unsigned len = 2 * (wToUseIFFT - wToUseFFT);

  std::vector<complex> container(3 * len);
  complex *x = &container.front(), *buffer = x + len, *y = buffer + len;
  {
    unsigned i = 0;
    for (unsigned iend = p1len; i != iend; ++i) {
      x[i].R = p1[i];
      x[i].I = 0;
    }
    for (unsigned iend = len; i != iend; ++i) {
      x[i].R = 0;
      x[i].I = 0;
    }

    i = 0;
    for (unsigned iend = p2len; i != iend; ++i) {
      y[i].R = p2[i];
      y[i].I = 0;
    }
    for (unsigned iend = len; i != iend; ++i) {
      y[i].R = 0;
      y[i].I = 0;
    }
  }

  fftConv(x, y, buffer, len, wToUseFFT, wToUseIFFT);
  for (unsigned i = 0, iend = p1len + p2len - 1; i != iend; ++i)
    rst[i] = buffer[i].R;
}
// FFT work space

void convEqualDelta(std::vector<losstype>::iterator xst,
                    std::vector<probtype>::iterator xpst, unsigned xlen,
                    std::vector<losstype>::iterator yst,
                    std::vector<probtype>::iterator ypst, unsigned ylen,
                    std::vector<losstype>::iterator rstst,
                    std::vector<probtype>::iterator rstpst,
                    wToUseFFT *w = NULL) {
  losstype d;
  *rstst = *xst + *yst;
  if (xlen > 1)
    d = *(xst + 1) - *xst;
  else if (ylen > 1)
    d = *(yst + 1) - *yst;
  else {
    if (xlen == 0 || ylen == 0) return;
    *rstpst = *xpst * *ypst;
    return;
  }
  //*rstst=*xst+*yst;
  std::vector<losstype>::iterator i = rstst + 1, end = rstst + xlen + ylen - 2;
  for (; i != end; ++i) *i = *(i - 1) + d;
  *i = *(xst + xlen - 1) + *(yst + ylen - 1);

  if (w == NULL) {
    for (unsigned i = 0; i != xlen; ++i) {
      for (unsigned j = 0; j != ylen; ++j) {
        rstpst[i + j] += xpst[i] * ypst[j];
      }
    }
    return;
  }

  fftConvFromReal(&*rstpst, &*xpst, &*ypst, xlen, ylen, *w);
}

void convEqualDelta(std::vector<probtype>::iterator xpst, unsigned xlen,
                    std::vector<probtype>::iterator ypst, unsigned ylen,
                    std::vector<probtype>::iterator rstpst,
                    wToUseFFT *w = NULL) {
  if (xlen == 0 || ylen == 0) return;
  if (w == NULL) {
    for (unsigned i = 0; i != xlen; ++i) {
      for (unsigned j = 0; j != ylen; ++j) {
        *(rstpst + i + j) += *(xpst + i) * *(ypst + j);
      }
    }
    return;
  }

  fftConvFromReal(&*rstpst, &*xpst, &*ypst, xlen, ylen, *w);
}

bool conv4productsCppAtomUntouched(
    std::vector<losstype>::iterator xst, std::vector<probtype>::iterator xpst,
    unsigned xlen, std::vector<losstype>::iterator yst,
    std::vector<probtype>::iterator ypst, unsigned ylen,
    std::vector<losstype>::iterator rstst,
    std::vector<probtype>::iterator rstpst, unsigned rstlen,
    unsigned forDeltaCompare, std::vector<losstype> &rst,
    std::vector<probtype> &rstp, rgMethod rglinear_add, wToUseFFT *w) {
  bool varianceEnlarged = 0;

  losstype /*xd=*(xst+1)-*xst, yd=*(yst+1)-*yst,*/
      xd = (*(xst + xlen - 1) - *xst) / (xlen - 1),
      yd = (*(yst + ylen - 1) - *yst) / (ylen - 1),
      rstd = ((*(xst + xlen - 1) - *xst) + (*(yst + ylen - 1) - *yst)) /
             (forDeltaCompare - 1);
  // std::cout<<"1.0----"<<forDeltaCompare<<"\n";
  losstype maxd = std::max(xd, std::max(yd, rstd));
  losstype ELIPSON = std::numeric_limits<losstype>::epsilon();
  if (std::abs(maxd - xd) < ELIPSON || std::abs(maxd - yd) < ELIPSON) {
    if (std::abs(maxd - yd) < ELIPSON) {
      std::swap(xst, yst);
      std::swap(xpst, ypst);
      std::swap(xlen, ylen);
    }  // now x is the one with greater delta. We need to regrid (y,yp) to x
    std::vector<losstype> ynew(lowerInt((*(yst + ylen - 2) - *yst) / maxd) + 3);
    std::vector<probtype> ypnew(ynew.size(), 0);
    formGrid(ynew.begin(), ynew.size() - 2, *yst, maxd);
    *(ynew.end() - 2) = *(yst + ylen - 2);
    ynew.back() = *(yst + ylen - 1);

    varianceEnlarged =
        rglinear_add(yst, ypst, ylen, ynew.begin(), ypnew.begin(), ynew.size());

    std::vector<losstype> conv(xlen + ynew.size());
    std::vector<probtype> convP(conv.size(), 0);  // doesn't need to be this
                                                  // long

    // unsigned finalrstlen;
    //--formGrid(rstst, *xst+*yst, *(yst+ylen)+*(xst+xlen), maxd, finalrstlen);
    // formGrid(rst.begin(), *xst+*yst, *(yst+ylen)+*(xst+xlen), maxd,
    // finalrstlen);

    convEqualDelta(xst, xpst, xlen, ynew.begin(), ypnew.begin(),
                   ynew.size() - 2, conv.begin(), convP.begin(), w);
    // rglinear_add(conv.begin(),convP.begin(),xlen,rst.begin(),rstp.begin(),finalrstlen);
    //--convEqualDelta(xst,xpst,xlen,ynew.begin(),ypnew.begin(),ynew.size()-2,rstst,rstpst);
    // convEqualDelta(xst,xpst,xlen,ynew.begin(),ypnew.begin(),ynew.size()-2,rst.begin(),rstp.begin());

    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), conv.size(),
                                    rstst, rstpst, rstlen);

    constMulV(xpst, xlen, convP.begin(), *(ypnew.end() - 2));
    constAddV(xst, xlen, conv.begin(), *(ynew.end() - 2));
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), xlen, rstst, rstpst, rstlen);
    // rglinear_add(conv.begin(),convP.begin(),xlen,rst.begin(),rstp.begin(),finalrstlen);
    constMulV(xpst, xlen, convP.begin(), *(ypnew.end() - 1));
    constAddV(xst, xlen, conv.begin(), *(ynew.end() - 1));
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), xlen, rstst, rstpst, rstlen);
    // rglinear_add(conv.begin(),convP.begin(),xlen,rst.begin(),rstp.begin(),finalrstlen);
    // rst.resize(finalrstlen);
    // rstp.resize(finalrstlen);
  } else {
    std::vector<losstype> xnew(
        lowerInt(((*(xst + xlen - 2) - *xst) / maxd) + 3)),
        ynew(lowerInt(((*(yst + ylen - 2) - *yst) / maxd) + 3));

    std::vector<probtype> xpnew(xnew.size(), 0), ypnew(ynew.size(), 0);

    formGrid(xnew.begin(), xnew.size() - 2, *xst, maxd);
    *(xnew.end() - 2) = *(xst + xlen - 2);
    xnew.back() = *(xst + xlen - 1);

    varianceEnlarged =
        rglinear_add(xst, xpst, xlen, xnew.begin(), xpnew.begin(), xnew.size());

    formGrid(ynew.begin(), ynew.size() - 2, *yst, maxd);
    *(ynew.end() - 2) = *(yst + ylen - 2);
    ynew.back() = *(yst + ylen - 1);

    varianceEnlarged =
        rglinear_add(yst, ypst, ylen, ynew.begin(), ypnew.begin(), ynew.size());

    std::vector<losstype> conv(xlen + ynew.size());
    std::vector<probtype> convP(conv.size(), 0);

    convEqualDelta(xnew.begin(), xpnew.begin(), xnew.size() - 2, ynew.begin(),
                   ypnew.begin(), ynew.size() - 2, conv.begin(), convP.begin(),
                   w);

    // unsigned finalrstlen;
    // formGrid(rst.begin(), *xst+*yst, *(yst+ylen)+*(xst+xlen), maxd,
    // finalrstlen);

    // convEqualDelta(xnew.begin(),xpnew.begin(),xnew.size()-2,ynew.begin(),ypnew.begin(),
    //              ynew.size()-2,rst.begin(),rstp.begin());

    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), conv.size(),
                                    rstst, rstpst, rstlen);

    constMulV(ypnew.begin(), ypnew.size(), convP.begin(), *(xpnew.end() - 2));
    constAddV(ynew.begin(), ynew.size(), conv.begin(), *(xnew.end() - 2));
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), ynew.size(),
                                    rstst, rstpst, rstlen);
    // rglinear_add(conv.begin(),convP.begin(),ynew.size(),rst.begin(),rstp.begin(),finalrstlen);

    constMulV(ypnew.begin(), ypnew.size(), convP.begin(), xpnew.back());
    constAddV(ynew.begin(), ynew.size(), conv.begin(), xnew.back());
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), ynew.size(),
                                    rstst, rstpst, rstlen);
    // rglinear_add(conv.begin(),convP.begin(),ynew.size(),rst.begin(),rstp.begin(),finalrstlen);

    constMulV(xpnew.begin(), xpnew.size() - 2, convP.begin(),
              *(ypnew.end() - 2));
    constAddV(xnew.begin(), xnew.size() - 2, conv.begin(), *(ynew.end() - 2));
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(),
                                    xnew.size() - 2, rstst, rstpst, rstlen);

    constMulV(xpnew.begin(), xpnew.size() - 2, convP.begin(), ypnew.back());
    constAddV(xnew.begin(), xnew.size() - 2, conv.begin(), ynew.back());
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(),
                                    xnew.size() - 2, rstst, rstpst, rstlen);
    // rst.resize(finalrstlen);
    // rstp.resize(finalrstlen);
  }
  return varianceEnlarged;
}

// [[Rcpp::export]]
DataFrame convSplitAtom4productsAtomUntouched(
    DataFrame X, DataFrame Y, String regridMethod = "lr", unsigned N = 256,
    unsigned forDeltaCompare = 256, double leftTrunc = 0, double rightTrunc = 0,
    bool useFFT = 0) {
  // std::vector<double> xR=X[0], xpR=X[1], yR=Y[0], ypR=Y[1];
  // std::vector<double> rst(N), rstp(N,0);

  std::vector<losstype> xR = X[0], yR = Y[0], rst(N);
  std::vector<probtype> xpR = X[1], ypR = Y[1], rstp(N, 0);

  truncateDist(xR, xpR, leftTrunc, rightTrunc);
  truncateDist(yR, ypR, leftTrunc, rightTrunc);

  formGrid(rst.begin(), rst.size() - 1, xR.front() + yR.front(),
           (xR.back() - xR.front() + yR.back() - yR.front()) / (N - 1));
  rst.back() = xR.back() + yR.back();

  if (xR.size() < 3 || yR.size() < 3) {
    std::cout << "One of the distributions has length less than 3";
    return Language("data.frame", 0, 0).eval();
  }

  rgMethod rglinear_add = nullptr;
  if (regridMethod == "lr") rglinear_add = &rglinearAddCpp;
  // else if(regridMethod=="r3")rglinear_add=&rglr3splitAddCpp;
  else if (regridMethod == "r4")
    rglinear_add = &rglr4splitAddCpp;
  else if (regridMethod == "lmm")
    rglinear_add = &lmmAddCpp;
  // else rglinear_add=&rglr4splitGaussianKernalAddCpp;
  else if (regridMethod == "lmm")
    rglinear_add = &lmmAddCpp;

  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);
  conv4productsCppAtomUntouched(xR.begin(), xpR.begin(), xR.size(), yR.begin(),
                                ypR.begin(), yR.size(), rst.begin(),
                                rstp.begin(), N, forDeltaCompare, rst, rstp,
                                rglinear_add, w);

  delete w;
  truncateDist(rst, rstp, leftTrunc, rightTrunc);
  return DataFrame::create(Named("val") = rst, Named("P") = rstp);
}

void printv(std::vector<double>::iterator st, unsigned len) {
  std::vector<double>::iterator end = st + len;
  for (; st != end; ++st) std::cout << *st << " ";
  std::cout << "\n";
}

bool conv9productsCppAtomUntouched(
    std::vector<losstype>::iterator xst, std::vector<probtype>::iterator xpst,
    unsigned xlen, std::vector<losstype>::iterator yst,
    std::vector<probtype>::iterator ypst, unsigned ylen,
    std::vector<losstype>::iterator rstst,
    std::vector<probtype>::iterator rstpst, unsigned rstlen,
    /*std::vector<double>&rst, std::vector<double>&rstp,*/
    unsigned forDeltaCompare, std::vector<losstype> &rst,
    std::vector<probtype> &rstp, rgMethod rglinear_add, wToUseFFT *w) {
  bool varianceEnlarged = 0;

  losstype xd = (*(xst + xlen - 1) - *xst) / (xlen - 1),
           yd = (*(yst + ylen - 1) - *yst) / (ylen - 1),
           rstd = ((*(xst + xlen - 1) - *xst) + (*(yst + ylen - 1) - *yst)) /
                  (forDeltaCompare - 1);
  // std::cout<<rstd<<"**\n";
  losstype maxd = std::max(xd, std::max(yd, rstd));
  losstype ELIPSON = std::numeric_limits<losstype>::epsilon();
  if (std::abs(maxd - xd) < ELIPSON || std::abs(maxd - yd) < ELIPSON) {
    if (std::abs(maxd - yd) < ELIPSON) {
      std::swap(xst, yst);
      std::swap(xpst, ypst);
      std::swap(xlen, ylen);
    }  // now x is the one with greater delta. We need to regrid (y,yp) to x
    std::vector<losstype> ynew(
        lowerInt((*(yst + ylen - 2) - *(yst + 1)) / maxd) + 4);
    std::vector<probtype> ypnew(ynew.size(), 0);

    ynew.front() = *yst;
    formGrid(ynew.begin() + 1, ynew.size() - 3, *(yst + 1), maxd);
    *(ynew.end() - 2) = *(yst + ylen - 2);
    ynew.back() = *(yst + ylen - 1);

    varianceEnlarged =
        rglinear_add(yst, ypst, ylen, ynew.begin(), ypnew.begin(), ynew.size());

    std::vector<losstype> conv(xlen + ynew.size());
    std::vector<probtype> convP(conv.size(), 0);

    convEqualDelta(xst, xpst, xlen, ynew.begin() + 1, ypnew.begin() + 1,
                   ynew.size() - 3, conv.begin(), convP.begin(), w);
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), conv.size(),
                                    rstst, rstpst, rstlen);

    constAddV(xst, xlen, conv.begin(), ynew.front());
    constMulV(xpst, xlen, convP.begin(), ypnew.front());
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), xlen, rstst, rstpst, rstlen);

    constAddV(xst, xlen, conv.begin(), *(ynew.end() - 2));
    constMulV(xpst, xlen, convP.begin(), *(ypnew.end() - 2));
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), xlen, rstst, rstpst, rstlen);

    constAddV(xst, xlen, conv.begin(), *(ynew.end() - 1));
    constMulV(xpst, xlen, convP.begin(), *(ypnew.end() - 1));
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), xlen, rstst, rstpst, rstlen);
  } else {
    std::vector<losstype> ynew(
        lowerInt((*(yst + ylen - 2) - *(yst + 1)) / maxd) + 4),
        xnew(lowerInt((*(xst + xlen - 2) - *(xst + 1)) / maxd) + 4);

    std::vector<probtype> ypnew(ynew.size(), 0), xpnew(xnew.size(), 0);

    ynew.front() = *yst;
    formGrid(ynew.begin() + 1, ynew.size() - 3, *(yst + 1), maxd);
    *(ynew.end() - 2) = *(yst + ylen - 2);
    ynew.back() = *(yst + ylen - 1);
    varianceEnlarged =
        rglinear_add(yst, ypst, ylen, ynew.begin(), ypnew.begin(), ynew.size());

    xnew.front() = *xst;
    formGrid(xnew.begin() + 1, xnew.size() - 3, *(xst + 1), maxd);
    *(xnew.end() - 2) = *(xst + xlen - 2);
    xnew.back() = *(xst + xlen - 1);

    varianceEnlarged =
        rglinear_add(xst, xpst, xlen, xnew.begin(), xpnew.begin(), xnew.size());

    std::vector<losstype> conv(xnew.size() + ynew.size());
    std::vector<probtype> convP(conv.size(), 0);

    convEqualDelta(xnew.begin() + 1, xpnew.begin() + 1, xnew.size() - 3,
                   ynew.begin() + 1, ypnew.begin() + 1, ynew.size() - 3,
                   conv.begin(), convP.begin(), w);

    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), xnew.size() + ynew.size() - 7,
                     rstst, rstpst, rstlen);

    constAddV(xnew.begin(), xnew.size(), conv.begin(), ynew.front());
    constMulV(xpnew.begin(), xpnew.size(), convP.begin(), ypnew.front());
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), xnew.size(),
                                    rstst, rstpst, rstlen);

    constAddV(xnew.begin(), xnew.size(), conv.begin(), *(ynew.end() - 2));
    constMulV(xpnew.begin(), xpnew.size(), convP.begin(), *(ypnew.end() - 2));
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), xnew.size(),
                                    rstst, rstpst, rstlen);

    constAddV(xnew.begin(), xnew.size(), conv.begin(), *(ynew.end() - 1));
    constMulV(xpnew.begin(), xpnew.size(), convP.begin(), *(ypnew.end() - 1));
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), xnew.size(),
                                    rstst, rstpst, rstlen);

    constAddV(ynew.begin() + 1, ynew.size() - 3, conv.begin(), xnew.front());
    constMulV(ypnew.begin() + 1, ypnew.size() - 3, convP.begin(),
              xpnew.front());
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(),
                                    ynew.size() - 3, rstst, rstpst, rstlen);

    constAddV(ynew.begin() + 1, ynew.size() - 3, conv.begin(),
              *(xnew.end() - 2));
    constMulV(ypnew.begin() + 1, ypnew.size() - 3, convP.begin(),
              *(xpnew.end() - 2));
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(),
                                    ynew.size() - 3, rstst, rstpst, rstlen);

    constAddV(ynew.begin() + 1, ynew.size() - 3, conv.begin(),
              *(xnew.end() - 1));
    constMulV(ypnew.begin() + 1, ypnew.size() - 3, convP.begin(),
              *(xpnew.end() - 1));
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(),
                                    ynew.size() - 3, rstst, rstpst, rstlen);
  }
  return varianceEnlarged;
}

// [[Rcpp::export]]
DataFrame convSplitAtom9productsAtomUntouched(
    DataFrame X, DataFrame Y, String regridMethod = "lr", unsigned N = 256,
    unsigned forDeltaCompare = 256, double leftTrunc = 0, double rightTrunc = 0,
    bool useFFT = 0) {
  // std::vector<double> xR=X[0], xpR=X[1], yR=Y[0], ypR=Y[1];
  // std::vector<double> rst(N), rstp(N);

  std::vector<losstype> xR = X[0], yR = Y[0], rst(N);
  std::vector<probtype> xpR = X[1], ypR = Y[1], rstp(N);

  truncateDist(xR, xpR, leftTrunc, rightTrunc);
  truncateDist(yR, ypR, leftTrunc, rightTrunc);

  formGrid(rst.begin(), rst.size() - 1, xR.front() + yR.front(),
           (xR.back() - xR.front() + yR.back() - yR.front()) / (N - 1));
  rst.back() = xR.back() + yR.back();

  if (xR.size() < 3 || yR.size() < 3) {
    std::cout << "One of the distributions has length less than 3";
    return Language("data.frame", 0, 0).eval();
  }
  // unsigned finalSize=N;
  rgMethod rglinear_add = nullptr;
  if (regridMethod == "lr") rglinear_add = &rglinearAddCpp;
  // else if(regridMethod=="r3")rglinear_add=&rglr3splitAddCpp;
  else if (regridMethod == "r4")
    rglinear_add = &rglr4splitAddCpp;
  else if (regridMethod == "lmm")
    rglinear_add = &lmmAddCpp;
  else if (regridMethod == "lmm")
    rglinear_add = &lmmAddCpp;
  // else rglinear_add=&rglr4splitGaussianKernalAddCpp;

  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);

  conv9productsCppAtomUntouched(xR.begin(), xpR.begin(), xR.size(), yR.begin(),
                                ypR.begin(), yR.size(), rst.begin(),
                                rstp.begin(), rst.size(), forDeltaCompare, rst,
                                rstp, rglinear_add, w);

  delete w;
  truncateDist(rst, rstp, leftTrunc, rightTrunc);
  return DataFrame::create(Named("val") = rst, Named("P") = rstp);
}

//----------------------------------------------------
bool conv4productsCppTouchAtom(
    std::vector<losstype>::iterator xst, std::vector<probtype>::iterator xpst,
    unsigned xlen, std::vector<losstype>::iterator yst,
    std::vector<probtype>::iterator ypst, unsigned ylen,
    std::vector<losstype>::iterator rstst,
    std::vector<probtype>::iterator rstpst, unsigned rstlen,
    unsigned forDeltaCompare, std::vector<losstype> &rst,
    std::vector<probtype> &rstp, rgMethod rglinear_add, wToUseFFT *w) {
  bool varianceEnlarged = 0;

  losstype xd = *(xst + 1) - *xst,
           yd = *(yst + 1) - *yst, /* xd=(*(xst+xlen-1)-*xst)/(xlen-1),
                                      yd=(*(yst+ylen-1)-*yst)/(ylen-1),*/
      rstd = ((*(xst + xlen - 1) - *xst) + (*(yst + ylen - 1) - *yst)) /
             (forDeltaCompare - 1);
  // std::cout<<"1.0----"<<forDeltaCompare<<"\n";
  losstype maxd = std::max(xd, std::max(yd, rstd));
  losstype ELIPSON = std::numeric_limits<losstype>::epsilon();
  if (std::abs(maxd - xd) < ELIPSON || std::abs(maxd - yd) < ELIPSON) {
    if (std::abs(maxd - yd) < ELIPSON) {
      std::swap(xst, yst);
      std::swap(xpst, ypst);
      std::swap(xlen, ylen);
    }  // now x is the one with greater delta. We need to regrid (y,yp) to x
    // std::cout<<"maxd=="<<maxd<<"\n";
    std::vector<losstype> ynew(lowerInt((*(yst + ylen - 1) - *yst) / maxd) + 2);
    std::vector<probtype> ypnew(ynew.size(), 0);
    // std::cout<<"y's size="<<ynew.size()<<"\n";
    formGrid(ynew.begin(), ynew.size() - 1, *yst, maxd);
    ynew.back() = *(yst + ylen - 1);

    varianceEnlarged =
        rglinear_add(yst, ypst, ylen, ynew.begin(), ypnew.begin(), ynew.size());

    std::vector<losstype> conv(xlen + ynew.size());
    std::vector<probtype> convP(conv.size(), 0);
    convEqualDelta(xst, xpst, xlen, ynew.begin(), ypnew.begin(),
                   ynew.size() - 1, conv.begin(), convP.begin(), w);

    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), conv.size(),
                                    rstst, rstpst, rstlen);

    constMulV(xpst, xlen, convP.begin(), *(ypnew.end() - 1));
    constAddV(xst, xlen, conv.begin(), *(ynew.end() - 1));

    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), xlen, rstst, rstpst, rstlen);
  } else {
    std::vector<losstype> xnew(
        lowerInt(((*(xst + xlen - 1) - *xst) / maxd) + 2)),
        ynew(lowerInt(((*(yst + ylen - 1) - *yst) / maxd) + 2));

    std::vector<probtype> xpnew(xnew.size(), 0), ypnew(ynew.size(), 0);

    formGrid(xnew.begin(), xnew.size() - 1, *xst, maxd);
    xnew.back() = *(xst + xlen - 1);
    varianceEnlarged =
        rglinear_add(xst, xpst, xlen, xnew.begin(), xpnew.begin(), xnew.size());

    formGrid(ynew.begin(), ynew.size() - 1, *yst, maxd);
    ynew.back() = *(yst + ylen - 1);
    varianceEnlarged =
        rglinear_add(yst, ypst, ylen, ynew.begin(), ypnew.begin(), ynew.size());

    std::vector<losstype> conv(xnew.size() + ynew.size(), 0);
    std::vector<probtype> convP(conv.size(), 0);

    convEqualDelta(xnew.begin(), xpnew.begin(), xnew.size() - 1, ynew.begin(),
                   ypnew.begin(), ynew.size() - 1, conv.begin(), convP.begin(),
                   w);
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), conv.size(),
                                    rstst, rstpst, rstlen);

    constMulV(ypnew.begin(), ypnew.size(), convP.begin(), *(xpnew.end() - 1));
    constAddV(ynew.begin(), ynew.size(), conv.begin(), *(xnew.end() - 1));
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), ynew.size(),
                                    rstst, rstpst, rstlen);

    constMulV(xpnew.begin(), xpnew.size() - 1, convP.begin(),
              *(ypnew.end() - 1));
    constAddV(xnew.begin(), xnew.size() - 1, conv.begin(), *(ynew.end() - 1));
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(),
                                    xnew.size() - 1, rstst, rstpst, rstlen);
  }
  return varianceEnlarged;
}

// [[Rcpp::export]]
DataFrame convSplitAtom4products(DataFrame X, DataFrame Y,
                                 String regridMethod = "lr", unsigned N = 256,
                                 unsigned forDeltaCompare = 256,
                                 double leftTrunc = 0, double rightTrunc = 0,
                                 bool useFFT = 0) {
  // std::vector<double> xR=X[0], xpR=X[1], yR=Y[0], ypR=Y[1];
  // std::vector<double> rst(N), rstp(N);

  std::vector<losstype> xR = X[0], yR = Y[0], rst(N);
  std::vector<probtype> xpR = X[1], ypR = Y[1], rstp(N);

  truncateDist(xR, xpR, leftTrunc, rightTrunc);
  truncateDist(yR, ypR, leftTrunc, rightTrunc);

  formGrid(rst.begin(), rst.size() - 1, xR.front() + yR.front(),
           (xR.back() - xR.front() + yR.back() - yR.front()) / (N - 1));
  rst.back() = xR.back() + yR.back();

  if (xR.size() < 3 || yR.size() < 3) {
    std::cout << "One of the distributions has length less than 3";
    return Language("data.frame", 0, 0).eval();
  }
  // unsigned finalSize=N;
  // std::cout<<forDeltaCompare<<"****\n";
  rgMethod rglinear_add = nullptr;
  if (regridMethod == "lr") rglinear_add = &rglinearAddCpp;
  // else if(regridMethod=="r3")rglinear_add=&rglr3splitAddCpp;
  else if (regridMethod == "r4")
    rglinear_add = &rglr4splitAddCpp;
  else if (regridMethod == "lmm")
    rglinear_add = &lmmAddCpp;
  // else rglinear_add=&rglr4splitGaussianKernalAddCpp;

  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);

  conv4productsCppTouchAtom(xR.begin(), xpR.begin(), xR.size(), yR.begin(),
                            ypR.begin(), yR.size(), rst.begin(), rstp.begin(),
                            rst.size(), forDeltaCompare, rst, rstp,
                            rglinear_add, w);
  delete w;
  truncateDist(rst, rstp, leftTrunc, rightTrunc);
  return DataFrame::create(Named("val") = rst, Named("P") = rstp);
}

bool conv9productsCppTouchAtom(
    std::vector<losstype>::iterator xst, std::vector<probtype>::iterator xpst,
    unsigned xlen, std::vector<losstype>::iterator yst,
    std::vector<probtype>::iterator ypst, unsigned ylen,
    std::vector<losstype>::iterator rstst,
    std::vector<probtype>::iterator rstpst, unsigned rstlen,
    unsigned forDeltaCompare, std::vector<losstype> &rst,
    std::vector<probtype> &rstp, rgMethod rglinear_add, wToUseFFT *w) {
  bool varianceEnlarged = 0;

  double xd = *(xst + 1) - *(xst),
         yd = *(yst + 1) - *(yst), /*xd=(*(xst+xlen-1)-*xst)/(xlen-1),
                                      yd=(*(yst+ylen-1)-*yst)/(ylen-1),*/
      rstd = ((*(xst + xlen - 1) - *xst) + (*(yst + ylen - 1) - *yst)) /
             (forDeltaCompare - 1);
  // std::cout<<rstd<<"**\n";
  double maxd = std::max(xd, std::max(yd, rstd));
  losstype ELIPSON = std::numeric_limits<losstype>::epsilon();
  if (std::abs(maxd - xd) < ELIPSON || std::abs(maxd - yd) < ELIPSON) {
    if (std::abs(maxd - yd) < ELIPSON) {
      std::swap(xst, yst);
      std::swap(xpst, ypst);
      std::swap(xlen, ylen);
    }  // now x is the one with greater delta. We need to regrid (y,yp) to x

    // std::vector<double>ynew(lowerInt((*(yst+ylen-1)-*(yst+1))/maxd)+3),
    // ypnew(ynew.size(),0);
    std::vector<losstype> ynew(
        lowerInt((*(yst + ylen - 1) - *(yst + 1)) / maxd) + 3);
    std::vector<probtype> ypnew(ynew.size(), 0);
    ynew.front() = *yst;
    formGrid(ynew.begin() + 1, ynew.size() - 2, *(yst + 1), maxd);
    ynew.back() = *(yst + ylen - 1);
    varianceEnlarged =
        rglinear_add(yst, ypst, ylen, ynew.begin(), ypnew.begin(), ynew.size());

    std::vector<losstype> conv(xlen + ynew.size());
    std::vector<probtype> convP(conv.size(), 0);
    convEqualDelta(xst, xpst, xlen, ynew.begin() + 1, ypnew.begin() + 1,
                   ynew.size() - 2, conv.begin(), convP.begin(), w);
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), conv.size(),
                                    rstst, rstpst, rstlen);

    constAddV(xst, xlen, conv.begin(), ynew.front());
    constMulV(xpst, xlen, convP.begin(), ypnew.front());
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), xlen, rstst, rstpst, rstlen);

    constAddV(xst, xlen, conv.begin(), *(ynew.end() - 1));
    constMulV(xpst, xlen, convP.begin(), *(ypnew.end() - 1));
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), xlen, rstst, rstpst, rstlen);
  } else {
    // std::cout<<1.1<<"\n";
    std::vector<losstype> ynew(
        lowerInt((*(yst + ylen - 1) - *(yst + 1)) / maxd) + 3);
    std::vector<probtype> ypnew(ynew.size(), 0);
    std::vector<losstype> xnew(
        lowerInt((*(xst + xlen - 1) - *(xst + 1)) / maxd) + 3);
    std::vector<probtype> xpnew(xnew.size(), 0);

    ynew.front() = *yst;
    formGrid(ynew.begin() + 1, ynew.size() - 2, *(yst + 1), maxd);
    ynew.back() = *(yst + ylen - 1);
    varianceEnlarged =
        rglinear_add(yst, ypst, ylen, ynew.begin(), ypnew.begin(), ynew.size());

    xnew.front() = *xst;
    formGrid(xnew.begin() + 1, xnew.size() - 2, *(xst + 1), maxd);
    xnew.back() = *(xst + xlen - 1);
    varianceEnlarged =
        rglinear_add(xst, xpst, xlen, xnew.begin(), xpnew.begin(), xnew.size());

    std::vector<losstype> conv(xnew.size() + ynew.size());
    std::vector<probtype> convP(conv.size(), 0);

    convEqualDelta(xnew.begin() + 1, xpnew.begin() + 1, xnew.size() - 2,
                   ynew.begin() + 1, ypnew.begin() + 1, ynew.size() - 2,
                   conv.begin(), convP.begin(), w);
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), xnew.size() + ynew.size() - 5,
                     rstst, rstpst, rstlen);

    constAddV(xnew.begin(), xnew.size(), conv.begin(), ynew.front());
    constMulV(xpnew.begin(), xpnew.size(), convP.begin(), ypnew.front());
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), xnew.size(),
                                    rstst, rstpst, rstlen);

    constAddV(xnew.begin(), xnew.size(), conv.begin(), *(ynew.end() - 1));
    constMulV(xpnew.begin(), xpnew.size(), convP.begin(), *(ypnew.end() - 1));
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), xnew.size(),
                                    rstst, rstpst, rstlen);

    constAddV(ynew.begin() + 1, ynew.size() - 2, conv.begin(), xnew.front());
    constMulV(ypnew.begin() + 1, ypnew.size() - 2, convP.begin(),
              xpnew.front());
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(),
                                    ynew.size() - 2, rstst, rstpst, rstlen);

    constAddV(ynew.begin() + 1, ynew.size() - 2, conv.begin(),
              *(xnew.end() - 1));
    constMulV(ypnew.begin() + 1, ypnew.size() - 2, convP.begin(),
              *(xpnew.end() - 1));
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(),
                                    ynew.size() - 2, rstst, rstpst, rstlen);
  }
  return varianceEnlarged;
}

// [[Rcpp::export]]
DataFrame convSplitAtom9products(DataFrame X, DataFrame Y,
                                 String regridMethod = "lr", unsigned N = 256,
                                 unsigned forDeltaCompare = 256,
                                 double leftTrunc = 0, double rightTrunc = 0,
                                 bool useFFT = 0) {
  std::vector<losstype> xR = X[0], yR = Y[0];
  std::vector<probtype> xpR = X[1], ypR = Y[1];
  std::vector<losstype> rst(N);
  std::vector<probtype> rstp(N);

  truncateDist(xR, xpR, leftTrunc, rightTrunc);
  truncateDist(yR, ypR, leftTrunc, rightTrunc);

  formGrid(rst.begin(), rst.size() - 1, xR.front() + yR.front(),
           (xR.back() - xR.front() + yR.back() - yR.front()) / (N - 1));
  rst.back() = xR.back() + yR.back();

  if (xR.size() < 3 || yR.size() < 3) {
    std::cout << "One of the distributions has length less than 3";
    return Language("data.frame", 0, 0).eval();
  }
  // unsigned finalSize=N;
  rgMethod rglinear_add = nullptr;
  if (regridMethod == "lr") rglinear_add = &rglinearAddCpp;
  // else if(regridMethod=="r3")rglinear_add=&rglr3splitAddCpp;
  else if (regridMethod == "r4")
    rglinear_add = &rglr4splitAddCpp;
  else if (regridMethod == "lmm")
    rglinear_add = &lmmAddCpp;
  // else rglinear_add=&rglr4splitGaussianKernalAddCpp;

  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);

  conv9productsCppTouchAtom(xR.begin(), xpR.begin(), xR.size(), yR.begin(),
                            ypR.begin(), yR.size(), rst.begin(), rstp.begin(),
                            rst.size(), forDeltaCompare, rst, rstp,
                            rglinear_add, w);
  delete w;

  truncateDist(rst, rstp, leftTrunc, rightTrunc);
  return DataFrame::create(Named("val") = rst, Named("P") = rstp);
}

// double meanOfSeq(std::vector<double>::iterator x,
// std::vector<double>::iterator xp, unsigned siz)

bool conv4productsIrregularCpp(
    std::vector<losstype>::iterator xst, std::vector<probtype>::iterator xpst,
    unsigned xlen, std::vector<losstype>::iterator yst,
    std::vector<probtype>::iterator ypst, unsigned ylen,
    std::vector<losstype>::iterator rstst,
    std::vector<probtype>::iterator rstpst, unsigned rstlen,
    unsigned forDeltaCompare, std::vector<losstype> &rst,
    std::vector<probtype> &rstp, rgMethod rglinear_add, wToUseFFT *w) {
  bool varianceEnlarged = 0;

  double xd = (*(xst + xlen - 2) - *xst) / (xlen - 2),
         yd = (*(yst + ylen - 2) - *yst) / (ylen - 2), rstd = 0;
  double xRange = *(xst + xlen - 1) - *xst, yRange = *(yst + ylen - 1) - *yst;
  double sumRange = xRange + yRange;

  for (unsigned tmpN = forDeltaCompare - 1; tmpN != 0; --tmpN) {
    rstd = sumRange / tmpN;
    if (lowerInt(sumRange / rstd) + 2 <= forDeltaCompare) break;
  }

  // std::cout<<"1.0----"<<forDeltaCompare<<"\n";
  double maxd = std::max(xd, std::max(yd, rstd));
  // std::cout<<"4product maxd=="<<maxd<<"\n";
  double ELIPSON = std::numeric_limits<double>::epsilon();
  if (std::abs(maxd - xd) < ELIPSON || std::abs(maxd - yd) < ELIPSON) {
    if (std::abs(maxd - yd) < ELIPSON) {
      std::swap(xst, yst);
      std::swap(xpst, ypst);
      std::swap(xlen, ylen);
    }  // now x is the one with greater delta. We need to regrid (y,yp) to x
    // std::cout<<"maxd=="<<maxd<<"\n";
    std::vector<losstype> ynew(lowerInt((*(yst + ylen - 1) - *yst) / maxd) + 2);
    std::vector<probtype> ypnew(ynew.size(), 0);
    // std::cout<<"y's size="<<ynew.size()<<"\n";
    formGrid(ynew.begin(), ynew.size() - 1, *yst, maxd);
    ynew.back() = *(yst + ylen - 1);

    varianceEnlarged =
        rglinear_add(yst, ypst, ylen, ynew.begin(), ypnew.begin(), ynew.size());

    // std::cout<<std::inner_product(ynew.begin(),ynew.end(),ypnew.begin(),0.0)<<"\n";

    unsigned finalrstlen =
        lowerInt((*(xst + xlen - 1) + ynew.back() - *xst - *yst) / maxd) + 2;
    if (finalrstlen > forDeltaCompare) finalrstlen = forDeltaCompare;

    rst.resize(finalrstlen);
    rstp.resize(finalrstlen);
    formGrid(rstst, finalrstlen - 1, *xst + *yst, maxd);
    rst.back() = *(xst + xlen - 1) + ynew.back();

    std::vector<losstype> conv(xlen + ynew.size());
    std::vector<probtype> convP(conv.size(), 0);
    // convEqualDelta(xpst, xlen, ypnew.begin(), ynew.size()-1, rstp.begin());
    convEqualDelta(xpst, xlen - 1, ypnew.begin(), ynew.size() - 1, rstp.begin(),
                   w);

    constMulV(xpst, xlen, convP.begin(), *(ypnew.end() - 1));
    constAddV(xst, xlen, conv.begin(), *(ynew.end() - 1));

    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), xlen,
                                    rst.begin(), rstp.begin(), rst.size());

    //***
    constMulV(ypnew.begin(), ynew.size() - 1, convP.begin(),
              *(xpst + xlen - 1));
    constAddV(ynew.begin(), ynew.size() - 1, conv.begin(), *(xst + xlen - 1));
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), ynew.size() - 1, rst.begin(),
                     rstp.begin(), rst.size());
    //***
  } else {
    // std::cout<<"1.2\n";
    std::vector<losstype> xnew(
        lowerInt(((*(xst + xlen - 1) - *xst) / maxd) + 2)),
        ynew(lowerInt(((*(yst + ylen - 1) - *yst) / maxd) + 2));

    std::vector<probtype> xpnew(xnew.size(), 0), ypnew(ynew.size(), 0);

    formGrid(xnew.begin(), xnew.size() - 1, *xst, maxd);
    xnew.back() = *(xst + xlen - 1);
    varianceEnlarged =
        rglinear_add(xst, xpst, xlen, xnew.begin(), xpnew.begin(), xnew.size());

    formGrid(ynew.begin(), ynew.size() - 1, *yst, maxd);
    ynew.back() = *(yst + ylen - 1);
    varianceEnlarged =
        rglinear_add(yst, ypst, ylen, ynew.begin(), ypnew.begin(), ynew.size());

    unsigned finalrstlen =
        lowerInt((xnew.back() + ynew.back() - xnew.front() - ynew.front()) /
                 maxd) +
        2;
    if (finalrstlen > forDeltaCompare) finalrstlen = forDeltaCompare;

    //   double xR=*(xst+xlen-1)-*xst, yR=ynew.back()-ynew.front();
    //   unsigned finalrstlen=unsigned(xR/maxd)+1+unsigned(yR/maxd)+1-1+1+
    //     unsigned((xR-unsigned(xR/maxd)*maxd+yR-unsigned(yR/maxd)*maxd)/maxd);

    rst.resize(finalrstlen);
    rstp.resize(finalrstlen);
    formGrid(rst.begin(), finalrstlen - 1, xnew.front() + ynew.front(), maxd);
    rst.back() = xnew.back() + ynew.back();

    std::vector<losstype> conv(xnew.size() + ynew.size(), 0);
    std::vector<probtype> convP(conv.size(), 0);

    //   convEqualDelta(xnew.begin(),xpnew.begin(),xnew.size()-1,ynew.begin(),ypnew.begin(),
    //                  ynew.size()-1,conv.begin(),convP.begin());
    convEqualDelta(xpnew.begin(), xnew.size() - 1, ypnew.begin(),
                   ynew.size() - 1, rstp.begin(), w);

    constMulV(ypnew.begin(), ypnew.size(), convP.begin(), *(xpnew.end() - 1));
    constAddV(ynew.begin(), ynew.size(), conv.begin(), *(xnew.end() - 1));
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), ynew.size(),
                                    rst.begin(), rstp.begin(), rst.size());

    constMulV(xpnew.begin(), xpnew.size() - 1, convP.begin(),
              *(ypnew.end() - 1));
    constAddV(xnew.begin(), xnew.size() - 1, conv.begin(), *(ynew.end() - 1));
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), xnew.size() - 1, rst.begin(),
                     rstp.begin(), rst.size());
  }
  return varianceEnlarged;
}








//'
//' Split-atom Four-products Convolution
//'
//' Convolve two PMFs using four-products convolution. Each PMF's support should 
//' be regular (have equal step sizes) except for the last step. \code{"Irregular"} 
//' refers to the fact that the size of the last step can be less.
//' 
//' @inheritParams rglr
//' 
//' @param Y Same as \code{X}.
//'
//' @param regridMethod  Regriding method in the routine:
//' 
//' \code{"lr"}: linear regrid. Default.
//' 
//' \code{"r4"}: four-point regrid.
//' 
//' \code{"r3"}: three-point regrid. For author's personal use.
//' 
//' \code{"r4g"}: four-point regrid with Gaussian kernel. For author's personal use.
//' 
//' Before convolution, one or both the PMFs need to be regrided such that they have
//' equal step size (main delta) except for the last steps (minor delta),
//' which can be smaller. See \code{forDeltaCompare} for how to determine the 
//' delta of the convolution PMF.
//' 
//' @param N  Maximum points allowed on the result PMF's support. Default 256.
//' 
//' 
//' @param forDeltaCompare For computing the main delta on the convolution grid.
//' Let \code{x} and \code{y} be the supports of \code{X} and \code{Y}, 
//' and \code{dx} and \code{dy} be the main deltas. The function first computes
//' the convolution grid's width \code{w = max(x) + max(y) - (min(x) + min(y))},
//' and then let \code{dc = max( dx, dy, w / forDeltaCompare )} be the main delta
//' on the convolution grid. The minor delta is thus \code{w \% dc}
//' In software release \code{forDeltaCompare == N}. 
//' Do not set \code{forDeltaCompare > N}.
//' 
//' 
//' @param leftTrunc  Nonnegative left tail truncation threshold. If nonzero, 
//' the left tail with probabilities \code{< leftTrunc} are truncated. Default 1e-10.
//' 
//' @param rightTrunc  Nonnegative right tail truncation threshold. If nonzero, 
//' the right tail with probabilities \code{< rightTrunc} are truncated. Default 1e-10.
//' 
//' @param useFFT Use Fast Fourier Transform. Default FALSE.
//'
//' @inherit rglr details
//'
//' @return A 2-column data frame as the result PMF. The PMF's support is 
//' regular except for the last step size which could be less.
//' 
//' @example inst/examples/convSplitAtom4productsIrregular.R
//'
// [[Rcpp::export]]
DataFrame convSplitAtom4productsIrregular(
    DataFrame X, 
    DataFrame Y, 
    String regridMethod = "lr", 
    int N = 256,
    int forDeltaCompare = 256,
    double leftTrunc = 1e-10, 
    double rightTrunc = 1e-10,
    bool useFFT = false
  ) 
{
  checkPMFintegrity(X);
  checkPMFintegrity(Y);
  
  
  std::vector<losstype> xR = X[0], yR = Y[0], rst(N);
  std::vector<probtype> xpR = X[1], ypR = Y[1], rstp(N, 0);

  
  truncateDist(xR, xpR, leftTrunc, rightTrunc);
  truncateDist(yR, ypR, leftTrunc, rightTrunc);

  // formGrid(rst.begin(),rst.size()-1,xR.front()+yR.front(),(xR.back()-xR.front()+yR.back()-yR.front())/(N-1));
  // rst.back()=xR.back()+yR.back();

  if (xR.size() < 3 || yR.size() < 3) {
    Rcout << "One of the distributions has length less than 3";
    return Language("data.frame", 0, 0).eval();
  }
  // unsigned finalSize=N;
  // std::cout<<forDeltaCompare<<"****\n";
  rgMethod rglinear_add = nullptr;
  if (regridMethod == "lr") rglinear_add = &rglinearAddCpp;
  // else if(regridMethod=="r3")rglinear_add=&rglr3splitAddCpp;
  else if (regridMethod == "r4")
    rglinear_add = &rglr4splitAddCpp;
  else if (regridMethod == "lmm")
    rglinear_add = &lmmAddCpp;
  // else rglinear_add=&rglr4splitGaussianKernalAddCpp;

  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);

  conv4productsIrregularCpp(xR.begin(), xpR.begin(), xR.size(), yR.begin(),
                            ypR.begin(), yR.size(), rst.begin(), rstp.begin(),
                            rst.size(), forDeltaCompare, rst, rstp,
                            rglinear_add, w);

  delete w;

  truncateDist(rst, rstp, leftTrunc, rightTrunc);
  // std::cout<<rst.size()<<"\n";
  // std::cout<<rstp.size()<<"\n";
  return DataFrame::create(Named("val") = rst, Named("P") = rstp);
}


bool conv9productsIrregularCpp(
    std::vector<losstype>::iterator xst, std::vector<probtype>::iterator xpst,
    unsigned xlen, std::vector<losstype>::iterator yst,
    std::vector<probtype>::iterator ypst, unsigned ylen,
    std::vector<losstype>::iterator rstst,
    std::vector<probtype>::iterator rstpst, unsigned rstlen,
    unsigned forDeltaCompare, std::vector<losstype> &rst,
    std::vector<probtype> &rstp, rgMethod rglinear_add, wToUseFFT *w) {
  bool varianceEnlarged = 0;

  double rstd = 0, xd = (*(xst + xlen - 2) - *(xst + 1)) / (xlen - 3),
               yd = (*(yst + ylen - 2) - *(yst + 1)) / (ylen - 3),
               xMidD = *(xst + xlen - 1) - *(xst + 1),
               yMidD = *(yst + ylen - 1) - *(yst + 1);

  double middleD = xMidD + yMidD;

  for (unsigned tmpN = forDeltaCompare - 1; tmpN != 0; --tmpN) {
    rstd = middleD / tmpN;
    //   if(unsigned(xMidD/rstd)+3-2+unsigned(yMidD/rstd)+3-2-1+1+
    //     unsigned((*(xst+1)+*(yst+1)-(*xst+*yst))/rstd)+1+
    //     unsigned((xMidD-unsigned(xMidD/rstd)*rstd+(yMidD-unsigned(yMidD/rstd)*rstd))/rstd)<=forDeltaCompare)break;

    if (lowerInt(middleD / rstd) + 2 + 1 +
            lowerInt((*(xst + 1) + *(yst + 1) - (*xst + *yst)) / rstd) <=
        forDeltaCompare)
      break;
  }

  double maxd = std::max(xd, std::max(yd, rstd));
  // std::cout<<"9product maxd=="<<maxd<<"\n";
  double ELIPSON = std::numeric_limits<double>::epsilon();
  if (std::abs(maxd - xd) < ELIPSON || std::abs(maxd - yd) < ELIPSON) {
    if (std::abs(maxd - yd) < ELIPSON) {
      std::swap(xst, yst);
      std::swap(xpst, ypst);
      std::swap(xlen, ylen);
    }  // now x is the one with greater delta. We need to regrid (y,yp) to x

    unsigned tmp;

    std::vector<losstype> ynew(
        lowerInt((*(yst + ylen - 1) - *(yst + 1)) / maxd) + 3);
    std::vector<probtype> ypnew(ynew.size(), 0);
    ynew.front() = *yst;
    formGrid(ynew.begin() + 1, ynew.size() - 2, *(yst + 1), maxd);
    ynew.back() = *(yst + ylen - 1);
    varianceEnlarged =
        rglinear_add(yst, ypst, ylen, ynew.begin(), ypnew.begin(), ynew.size());

    unsigned finalrstlen;
    tmp = lowerInt((*(ynew.begin() + 1) + *(xst + 1) - (ynew.front() + *xst)) /
                   maxd);

    finalrstlen =
        lowerInt((*(xst + xlen - 1) + ynew.back() - (*(xst + 1) + *(yst + 1))) /
                 maxd) +
        3 + tmp;

    rst.resize(finalrstlen);
    rstp.resize(finalrstlen);

    rst.front() = *xst + ynew.front();
    //   for(std::vector<double>::iterator
    //   i=rst.begin()+1,end=rst.begin()+tmp+1;i!=end;++i)
    //     *i=*(i-1)+maxd;
    //   formGrid(rst.begin()+tmp+1, rst.end()-2-(rst.begin()+tmp+1)+1,
    //   *(xst+1)+*(ynew.begin()+1), maxd);
    double tmpNum = *(xst + 1) + *(ynew.begin() + 1) - tmp * maxd;
    if (tmpNum > 0)
      formGrid(rst.begin() + 1, finalrstlen - 2, tmpNum, maxd);
    else
      formGrid(rst.begin() + 1, finalrstlen - 2, rst.front() + maxd, maxd);

    rst.back() = *(xst + xlen - 1) + ynew.back();

    std::vector<losstype> conv(xlen + ynew.size());
    std::vector<probtype> convP(conv.size(), 0);

    // convEqualDelta(xst,xpst,xlen,ynew.begin()+1,ypnew.begin()+1,ynew.size()-2,conv.begin(),convP.begin());

    convEqualDelta(xpst + 1, xlen - 2, ypnew.begin() + 1, ynew.size() - 2,
                   rstp.begin() + tmp + 1, w);

    // varianceEnlarged=
    // rglinear_add(conv.begin(),convP.begin(),conv.size(),rstst,rstpst,rstlen);

    constAddV(xst, xlen, conv.begin(), ynew.front());
    constMulV(xpst, xlen, convP.begin(), ypnew.front());
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), xlen,
                                    rst.begin(), rstp.begin(), finalrstlen);

    constAddV(xst, xlen, conv.begin(), *(ynew.end() - 1));
    constMulV(xpst, xlen, convP.begin(), *(ypnew.end() - 1));
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), xlen,
                                    rst.begin(), rstp.begin(), finalrstlen);

    constAddV(ynew.begin() + 1, ynew.size() - 2, conv.begin(), *xst);
    constMulV(ypnew.begin() + 1, ypnew.size() - 2, convP.begin(), *xpst);
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), ynew.size() - 2, rst.begin(),
                     rstp.begin(), finalrstlen);

    constAddV(ynew.begin() + 1, ynew.size() - 2, conv.begin(),
              *(xst + xlen - 1));
    constMulV(ypnew.begin() + 1, ypnew.size() - 2, convP.begin(),
              *(xpst + xlen - 1));
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), ynew.size() - 2, rst.begin(),
                     rstp.begin(), finalrstlen);
  } 
  else 
  {
    std::vector<losstype> ynew(
        lowerInt((*(yst + ylen - 1) - *(yst + 1)) / maxd) + 3),
        xnew(lowerInt((*(xst + xlen - 1) - *(xst + 1)) / maxd) + 3);
    std::vector<probtype> ypnew(ynew.size(), 0), xpnew(xnew.size(), 0);

    ynew.front() = *yst;
    formGrid(ynew.begin() + 1, ynew.size() - 2, *(yst + 1), maxd);
    ynew.back() = *(yst + ylen - 1);
    varianceEnlarged =
        rglinear_add(yst, ypst, ylen, ynew.begin(), ypnew.begin(), ynew.size());

    xnew.front() = *xst;
    formGrid(xnew.begin() + 1, xnew.size() - 2, *(xst + 1), maxd);
    xnew.back() = *(xst + xlen - 1);
    varianceEnlarged =
        rglinear_add(xst, xpst, xlen, xnew.begin(), xpnew.begin(), xnew.size());

    unsigned finalrstlen;
    unsigned tmp = lowerInt(
        (*(ynew.begin() + 1) + *(xst + 1) - (ynew.front() + *xst)) / maxd);
    
    
    finalrstlen = lowerInt((xnew.back() + ynew.back() -
                            (*(xnew.begin() + 1) + *(ynew.begin() + 1))) /
                           maxd) + 3 + tmp;
    if (finalrstlen > forDeltaCompare) finalrstlen = forDeltaCompare;

    rst.resize(finalrstlen);
    rstp.resize(finalrstlen);

    rst.front() = xnew.front() + ynew.front();
    

    double tmpNum = *(xst + 1) + *(ynew.begin() + 1) - tmp * maxd;
    if (tmpNum > 0)
      formGrid(rst.begin() + 1, finalrstlen - 2, tmpNum, maxd);
    else
      formGrid(rst.begin() + 1, finalrstlen - 2, rst.front() + maxd, maxd);

    rst.back() = xnew.back() + ynew.back();

    std::vector<losstype> conv(xnew.size() + ynew.size());
    std::vector<probtype> convP(conv.size(), 0);

    
    convEqualDelta(xpnew.begin() + 1, xnew.size() - 2, ypnew.begin() + 1,
                   ynew.size() - 2, rstp.begin() + tmp + 1, w);

    
    constAddV(xnew.begin(), xnew.size(), conv.begin(), ynew.front());
    constMulV(xpnew.begin(), xpnew.size(), convP.begin(), ypnew.front());
    
    
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), xnew.size(),
                                    rst.begin(), rstp.begin(), finalrstlen);

    
    constAddV(xnew.begin(), xnew.size(), conv.begin(), *(ynew.end() - 1));
    constMulV(xpnew.begin(), xpnew.size(), convP.begin(), *(ypnew.end() - 1));
    
    
    varianceEnlarged = rglinear_add(conv.begin(), convP.begin(), xnew.size(),
                                    rst.begin(), rstp.begin(), finalrstlen);

    
    constAddV(ynew.begin() + 1, ynew.size() - 2, conv.begin(), xnew.front());
    constMulV(ypnew.begin() + 1, ypnew.size() - 2, convP.begin(),
              xpnew.front());
    
    
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), ynew.size() - 2, rst.begin(),
                     rstp.begin(), finalrstlen);

    constAddV(ynew.begin() + 1, ynew.size() - 2, conv.begin(),
              *(xnew.end() - 1));
    constMulV(ypnew.begin() + 1, ypnew.size() - 2, convP.begin(),
              *(xpnew.end() - 1));
    
    
    varianceEnlarged =
        rglinear_add(conv.begin(), convP.begin(), ynew.size() - 2, rst.begin(),
                     rstp.begin(), finalrstlen);
  }
  return varianceEnlarged;
}


// [[Rcpp::export]]
DataFrame convSplitAtom9productsIrregular(
    DataFrame X, DataFrame Y, String regridMethod = "lr", unsigned N = 256,
    unsigned forDeltaCompare = 256, double leftTrunc = 0, double rightTrunc = 0,
    bool useFFT = 0) 
{
  std::vector<losstype> xR = X[0], yR = Y[0], rst(N);
  std::vector<probtype> xpR = X[1], ypR = Y[1], rstp(N);

  truncateDist(xR, xpR, leftTrunc, rightTrunc);
  truncateDist(yR, ypR, leftTrunc, rightTrunc);

  
  if (xR.size() < 3 || yR.size() < 3) {
    std::cout << "One of the distributions has length less than 3";
    return Language("data.frame", 0, 0).eval();
  }
  rgMethod rglinear_add = nullptr;
  if (regridMethod == "lr") rglinear_add = &rglinearAddCpp;
  else if (regridMethod == "r4")
    rglinear_add = &rglr4splitAddCpp;
  else if (regridMethod == "lmm")
    rglinear_add = &lmmAddCpp;

  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);

  conv9productsIrregularCpp(xR.begin(), xpR.begin(), xR.size(), yR.begin(),
                            ypR.begin(), yR.size(), rst.begin(), rstp.begin(),
                            rst.size(), forDeltaCompare, rst, rstp,
                            rglinear_add, w);

  delete w;
  truncateDist(rst, rstp, leftTrunc, rightTrunc);
  return DataFrame::create(Named("val") = rst, Named("P") = rstp);
}


//***************************************************************************************************

// Test the speed and accuracy and everything..
bool convBtRgdForTest(std::vector<losstype>::iterator xst,
                      std::vector<probtype>::iterator xpst, unsigned xlen,
                      std::vector<losstype>::iterator yst,
                      std::vector<probtype>::iterator ypst, unsigned ylen,
                      std::vector<losstype>::iterator rstst,
                      std::vector<probtype>::iterator rstpst, unsigned rstlen,
                      unsigned forDeltaCompare, std::vector<losstype> &rst,
                      std::vector<probtype> &rstp, rgMethod rglinear_add,
                      wToUseFFT *w) {
  bool varianceEnlarged = 0;
  if (xlen > ylen) 
  {
    std::swap(xlen, ylen);
    std::swap(xst, yst);
    std::swap(xpst, ypst);
  }

  
  std::vector<losstype> ytmp(ylen);
  std::vector<probtype> yptmp(ylen);
  for (unsigned i = 0; i != xlen; ++i) {
    constAddV(yst, ylen, ytmp.begin(), *(xst + i));
    constMulV(ypst, ylen, yptmp.begin(), *(xpst + i));
    varianceEnlarged =
        rglinear_add(ytmp.begin(), yptmp.begin(), ylen, rstst, rstpst, rstlen);
  }
  return varianceEnlarged;
}


struct dist 
{
  std::vector<losstype> val;
  std::vector<probtype> P;
  dist(unsigned N) {
    val.resize(N);
    P.resize(N, 0);
  }
  dist(){};
  void getGrid(losstype start, losstype end) {
    val.front() = start;
    losstype d = (end - start) / (val.size() - 1);
    for (std::vector<losstype>::iterator i = val.begin() + 1; i != val.end();
         ++i)
      *i = *(i - 1) + d;
    val.back() = end;
    std::fill(P.begin(), P.end(), 0);
  }
};


void swapDist(dist &x, dist &y) 
{
  std::swap(x.val, y.val);
  std::swap(x.P, y.P);
}


struct pointerToDist 
{
  losstype range;
  std::vector<dist>::iterator ptr;
  bool operator<(const pointerToDist &str) const { return range < str.range; }
  bool operator<(const losstype rg) const { return range < rg; }
};


typedef bool (*convMethod)(std::vector<losstype>::iterator,
                           std::vector<probtype>::iterator, unsigned,
                           std::vector<losstype>::iterator,
                           std::vector<probtype>::iterator, unsigned,
                           std::vector<losstype>::iterator,
                           std::vector<probtype>::iterator, unsigned, unsigned,
                           std::vector<losstype> &, std::vector<probtype> &,
                           rgMethod, wToUseFFT *);


typedef bool (*truncMethod)(std::vector<losstype> &, std::vector<probtype> &,
                            double, double);


void sequentialConv(std::vector<pointerToDist> &order,
                    std::vector<losstype> &rst, std::vector<probtype> &rstp,
                    String &convolutionMethod, String &regridMethod,
                    double headTrunc, double tailTrunc, bool single, unsigned N,
                    unsigned forDeltaCompare, wToUseFFT *w) {
  convMethod conv;
  
  
  if (convolutionMethod == "splitAtom4products")
    conv = &conv4productsCppTouchAtom;
  else if (convolutionMethod == "splitAtom9products")
    conv = &conv9productsCppTouchAtom;
  else if (convolutionMethod == "splitAtom4productsAtomUntouched")
    conv = &conv4productsCppAtomUntouched;
  else if (convolutionMethod == "splitAtom9productsAtomUntouched")
    conv = &conv9productsCppAtomUntouched;
  else if (convolutionMethod == "splitAtom4productsIrregular")
    conv = &conv4productsIrregularCpp;
  else if (convolutionMethod == "splitAtom9productsIrregular")
    conv = &conv9productsIrregularCpp;
  else conv = &convBtRgdForTest;

  
  rgMethod rglinear_add = nullptr;
  if (regridMethod == "lr")       rglinear_add = &rglinearAddCpp;
  else if (regridMethod == "r4")  rglinear_add = &rglr4splitAddCpp;
  else if (regridMethod == "lmm") rglinear_add = &lmmAddCpp;
  

  truncMethod truncate;
  if (single) truncate = &truncateDistSingle;
  else truncate = &truncateDist;

  
  std::vector<dist>::iterator i = (order.begin() + 1)->ptr;

  
  formGrid(rst.begin(), rst.size(),
           order.begin()->ptr->val.front() + i->val.front(),
           (order.begin()->ptr->val.back() - order.begin()->ptr->val.front() +
            (i->val.back() - i->val.front())) /
               (rst.size() - 1));

  rst.back() = order.begin()->ptr->val.back() + i->val.back();

  
  conv(order.begin()->ptr->val.begin(), order.begin()->ptr->P.begin(),
       order.begin()->ptr->val.size(), i->val.begin(), i->P.begin(),
       i->val.size(), rst.begin(), rstp.begin(), rst.size(), forDeltaCompare,
       rst, rstp, rglinear_add, w);

  
  truncate(rst, rstp, headTrunc, tailTrunc);

  
  std::vector<losstype> accRstVal(N);
  std::vector<probtype> accRstP(N);

  
  std::vector<pointerToDist>::iterator k = order.begin() + 2;

  
  while (k != order.end()) {
    std::swap(rst, accRstVal);
    std::swap(rstp, accRstP);
    rst.resize(N);
    rstp.resize(N);
    formGrid(rst.begin(), rst.size(), accRstVal.front() + k->ptr->val.front(),
             (accRstVal.back() - accRstVal.front() +
              (k->ptr->val.back() - k->ptr->val.front())) /
                 (rst.size() - 1));
    rst.back() = accRstVal.back() + k->ptr->val.back();
    std::fill(rstp.begin(), rstp.end(), 0);
    conv(accRstVal.begin(), accRstP.begin(), accRstP.size(),
         k->ptr->val.begin(), k->ptr->P.begin(), k->ptr->val.size(),
         rst.begin(), rstp.begin(), rst.size(), forDeltaCompare, rst, rstp,
         rglinear_add, w);
    truncate(rst, rstp, headTrunc, tailTrunc);
    //  ++i;
    ++k;
  }
}




//' 
//' Sequential convolution
//' 
//' Convolve a list of PMFs sequentially.
//' 
//' @param lisOfDists 
//'   A list of PMFs. 
//' 
//'   Each PMF is a 2-column data frame \code{(support, probability)}. 
//'   \code{support} should be regular (equal step sizes). \code{probability} 
//'   should be nonnegative. 
//' 
//'   A degenerate PMF must be a 1-row data frame. A degenerate PMF represented as 
//'   \code{data.frame(val = c(0, 1, 2), P = c(0, 0, 1))} instead of 
//'   \code{data.frame(val = 2, P = 1)} can lead to runtime error.
//' 
//' 
//' @param convMethod
//'   Convolution method:
//'   
//'   \code{"splitAtom4productsIrregular"}: use 
//'   \code{\link{convSplitAtom4productsIrregular}}, the routine in 
//'   software release. Default.
//' 
//'   \code{"convbt"}: use \code{\link{convBt}}, brute-force convolution.
//' 
//'   The following were implemented for personal use:
//' 
//'   \code{"splitAtom4products"}: use \code{\link{convSplitAtom4products}}.
//' 
//'   \code{"splitAtom9products"}: use \code{\link{convSplitAtom9products}}.
//' 
//'   \code{"splitAtom4productsAtomUntouched"}: 
//'     use \code{\link{convSplitAtom4productsAtomUntouched}}.
//' 
//'   \code{"splitAtom9productsAtomUntouched"}: 
//'     use \code{\link{convSplitAtom9productsAtomUntouched}}.
//' 
//'   \code{"splitAtom9productsIrregular"}: 
//'     use \code{\link{convSplitAtom9productsIrregular}}.
//' 
//'
//' @inheritParams convSplitAtom4productsIrregular
//' 
//' @inherit convSplitAtom4productsIrregular return
//' 
//' @details See convolution function, e.g., \code{\link{convSplitAtom4productsIrregular}}
//' and \href{https://www.mdpi.com/2227-9091/7/2/54/htm}{Direct and Hierarchical 
//' Models for Aggregating Spatially Dependent Catastrophe Risks}.
//' 
//' @example inst/examples/sequentialConvAll.R
// [[Rcpp::export]]
DataFrame sequentialConvAll(List lisOfDists, 
                            String convMethod = "splitAtom4productsIrregular", 
                            double headTrunc = 1e-10, 
                            double tailTrunc = 1e-10,
                            String regridMethod = "lr", 
                            int N = 256,
                            int forDeltaCompare = 256, 
                            bool useFFT = false) 
{
  
  
  checkPMFlistIntegrity(lisOfDists);
  
  
  bool single = true;
  if (lisOfDists.size() == 1) 
  {
    List rst = lisOfDists[0];
    return rst;
  }

  std::vector<dist> lisOfDistsCpp(lisOfDists.size());
  int k = 0;
  double degenerateValSum = 0;
  for (int i = 0; i != lisOfDists.size(); ++i) {
    List tmp = lisOfDists[i];
    NumericVector value = tmp[0], prob = tmp[1];

    if (value.size() > 2) {
      lisOfDistsCpp[k].val.assign(value.begin(), value.end());
      lisOfDistsCpp[k].P.assign(prob.begin(), prob.end());
      k += 1;
    } 
    else if ( value.size() == 2 )
    {
      lisOfDistsCpp[k].val.resize(3);
      lisOfDistsCpp[k].val[0] = value[0];
      lisOfDistsCpp[k].val[1] = (value[0] + value[1]) / 2;
      lisOfDistsCpp[k].val[2] = value[1];
      lisOfDistsCpp[k].P.resize(3);
      lisOfDistsCpp[k].P[0] = prob[0];
      lisOfDistsCpp[k].P[1] = 0;
      lisOfDistsCpp[k].P[2] = prob[1];
      k += 1;
    }
    else degenerateValSum += value[0];
  }
  lisOfDistsCpp.resize(k);
  

  std::vector<losstype> rst(N);
  std::vector<probtype> rstp(N, 0);
  std::vector<pointerToDist> order(lisOfDistsCpp.size());
  std::vector<dist>::iterator i = lisOfDistsCpp.begin();
  std::vector<pointerToDist>::iterator j = order.begin();
  for (; i != lisOfDistsCpp.end(); ++i) {
    if (i->val.size() > 1) {
      j->range = (i->val.back() - i->val.front()) / (i->val.size() + 0.0);
      j->ptr = i;
      ++j;
    }
  }
  order.resize(j - order.begin());

  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);

  sequentialConv(order, rst, rstp, convMethod, regridMethod, headTrunc,
                 tailTrunc, single, N, forDeltaCompare, w);
  if (useFFT) delete w;

  double sumP = std::accumulate(rstp.begin(), rstp.end(), 0.0);
  for (std::vector<probtype>::iterator i = rstp.begin(); i != rstp.end(); ++i)
    *i /= sumP;
  
  
  NumericVector finalRst(rst.size());
  for (int i = 0, iend = rst.size(); i < iend; ++i)
    finalRst[i] = rst[i] + degenerateValSum;
  

  return DataFrame::create(Named("val") = finalRst, Named("P") = rstp);
}




//'
//' Sorted sequential convolution
//' 
//' Sort the list of distributions by distribution's support step size, and then 
//' call \code{\link{sequentialConvAll}}. This MAY mitigate variance inflation.
//' 
//' @inherit sequentialConvAll
//'
//' @example inst/examples/sortedSequentialConvAll.R
// [[Rcpp::export]]
DataFrame sortedSequentialConvAll(
    List lisOfDists,
    String convMethod = "splitAtom4productsIrregular",
    double headTrunc = 1e-10,
    double tailTrunc = 1e-10, 
    String regridMethod = "lr",
    int N = 256,
    int forDeltaCompare = 256,
    bool useFFT = false) 
{
  
  
  checkPMFlistIntegrity(lisOfDists);
  
  
  bool single = 1; // Truncation method.
  if (lisOfDists.size() == 1) {
    List rst = lisOfDists[0];
    return rst;
  }

  
  std::vector<dist> lisOfDistsCpp(lisOfDists.size());
  int k = 0;
  double degenerateValSum = 0;
  for (int i = 0; i != lisOfDists.size(); ++i) {
    List tmp = lisOfDists[i];
    NumericVector value = tmp[0], prob = tmp[1];

    if (value.size() > 2) 
    {
      lisOfDistsCpp[k].val.assign(value.begin(), value.end());
      lisOfDistsCpp[k].P.assign(prob.begin(), prob.end());
      k += 1;
    } 
    else if (value.size() == 2)
    {
      lisOfDistsCpp[k].val.resize(3);
      lisOfDistsCpp[k].val[0] = value[0];
      lisOfDistsCpp[k].val[1] = (value[0] + value[1]) / 2;
      lisOfDistsCpp[k].val[2] = value[1];
      lisOfDistsCpp[k].P.resize(3);
      lisOfDistsCpp[k].P[0] = prob[0];
      lisOfDistsCpp[k].P[1] = 0;
      lisOfDistsCpp[k].P[2] = prob[1];
      k += 1;
    }
    else degenerateValSum += value[0];
  }
  lisOfDistsCpp.resize(k);
  

  std::vector<losstype> rst(N);
  std::vector<probtype> rstp(N, 0);
  std::vector<pointerToDist> order(lisOfDistsCpp.size());
  std::vector<dist>::iterator i = lisOfDistsCpp.begin();
  std::vector<pointerToDist>::iterator j = order.begin();

  
  for (; i != lisOfDistsCpp.end(); ++i) 
  {
    if (i->val.size() > 1) 
    {
      j->range = i->val.back() - i->val.front();
      j->ptr = i;
      ++j;
    }
  }
  order.resize(j - order.begin());
  std::sort(order.begin(), order.end());

  
  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);
  sequentialConv(order, rst, rstp, convMethod, regridMethod, headTrunc,
                 tailTrunc, single, N, forDeltaCompare, w);
  if (useFFT) delete w;
  double sumP = std::accumulate(rstp.begin(), rstp.end(), 0.0);
  for (std::vector<probtype>::iterator i = rstp.begin(); i != rstp.end(); ++i)
    *i /= sumP;
  
  
  NumericVector finalRst(rst.size());
  for (int i = 0, iend = rst.size(); i < iend; ++i)
    finalRst[i] = rst[i] + degenerateValSum;
  
  
  return DataFrame::create(Named("val") = finalRst, Named("P") = rstp);
}


// first and last are pointint to orderAcc
void divide2sums(std::vector<pointerToDist>::iterator first,
                          std::vector<pointerToDist>::iterator last,
                          std::vector<pointerToDist>::iterator &head,
                          convMethod &conv, rgMethod &rglinear_add, unsigned &N,
                          unsigned &forDeltaCompare,
                          truncMethod &truncate, double headTrunc,
                          double tailTrunc, dist &functionValue, wToUseFFT *w) 
{
  if (first + 1 == last) 
  {
    dist rst(N);
    formGrid(rst.val.begin(), rst.val.size(),
             first->ptr->val.front() + last->ptr->val.front(),
             (first->ptr->val.back() - first->ptr->val.front() +
              (last->ptr->val.back() - last->ptr->val.front())) /
                 (rst.val.size() - 1));
    rst.val.back() = first->ptr->val.back() + last->ptr->val.back();
    conv(first->ptr->val.begin(), first->ptr->P.begin(), first->ptr->P.size(),
         last->ptr->val.begin(), last->ptr->P.begin(), last->ptr->P.size(),
         rst.val.begin(), rst.P.begin(), N, forDeltaCompare, rst.val, rst.P,
         rglinear_add, w);
    truncate(rst.val, rst.P, headTrunc, tailTrunc);
    swapDist(rst, functionValue);
  }
  

  pointerToDist sum;
  sum.range = (first->range + last->range) / 2;
  std::vector<pointerToDist>::iterator splitPoint =
      std::lower_bound(first, last + 1, sum);
  
  
  if (splitPoint->range == sum.range && splitPoint != last) {
    if (splitPoint->range - (splitPoint - 1)->range >
        (splitPoint + 1)->range - splitPoint->range)
      ++splitPoint;
  }

  
  dist leftBranch, rightBranch, &left = *(first->ptr), &right = *(last->ptr);

  
  if (first != splitPoint - 1) {
    divide2sums(first, splitPoint - 1, head, conv, rglinear_add,
                               N, forDeltaCompare, truncate,
                               headTrunc, tailTrunc, leftBranch, w);
    left = leftBranch;
  }
  if (splitPoint != last) {
    divide2sums(splitPoint, last, head, conv, rglinear_add, N,
                forDeltaCompare, truncate,
                headTrunc, tailTrunc, rightBranch, w);
    right = rightBranch;
  }

  dist rst(N);

  formGrid(rst.val.begin(), rst.val.size(),
           left.val.front() + right.val.front(),
           (left.val.back() - left.val.front() +
            (right.val.back() - right.val.front())) /
               (N - 1));
  rst.val.back() = left.val.back() + right.val.back();
  conv(left.val.begin(), left.P.begin(), left.P.size(), right.val.begin(),
       right.P.begin(), right.P.size(), rst.val.begin(), rst.P.begin(), N,
       forDeltaCompare, rst.val, rst.P, rglinear_add, w);
  truncate(rst.val, rst.P, headTrunc, tailTrunc);
  swapDist(rst, functionValue);
}


void closestPairConv(std::vector<pointerToDist> &order,
                     String &convolutionMethod, String &regridMethod,
                     unsigned N, unsigned forDeltaCompare, bool single,
                     double headTrunc, double tailTrunc,
                     dist &functionValue, wToUseFFT *w) {
  
  
  convMethod conv;
  
  
  if (convolutionMethod == "splitAtom4products")
    conv = &conv4productsCppTouchAtom;
  else if (convolutionMethod == "splitAtom9products")
    conv = &conv9productsCppTouchAtom;
  else if (convolutionMethod == "splitAtom4productsAtomUntouched")
    conv = &conv4productsCppAtomUntouched;
  else if (convolutionMethod == "splitAtom9productsAtomUntouched")
    conv = &conv9productsCppAtomUntouched;
  else if (convolutionMethod == "splitAtom4productsIrregular")
    conv = &conv4productsIrregularCpp;
  else if (convolutionMethod == "splitAtom9productsIrregular")
    conv = &conv9productsIrregularCpp;
  else conv = &convBtRgdForTest;

  
  rgMethod rglinear_add = nullptr;
  if (regridMethod == "lr") rglinear_add = &rglinearAddCpp;
  else if (regridMethod == "r4")
    rglinear_add = &rglr4splitAddCpp;
  else if (regridMethod == "lmm")
    rglinear_add = &lmmAddCpp;
  

  truncMethod truncate;
  if (single)
    truncate = &truncateDistSingle;
  else
    truncate = &truncateDist;

  std::vector<pointerToDist> orderAcc(order.begin(), order.end());
  for (std::vector<pointerToDist>::iterator i = orderAcc.begin() + 1;
       i != orderAcc.end(); ++i)
    i->range += (i - 1)->range;
  std::vector<pointerToDist>::iterator thehead = orderAcc.begin(),
                                       &head = thehead;
  divide2sums(orderAcc.begin(), orderAcc.end() - 1, head, conv,
              rglinear_add, N, forDeltaCompare, truncate,
              headTrunc, tailTrunc, functionValue, w);
}




//'
//' Closest pair convolution
//' 
//' Recursively searches for pairs of PMFs that have the closest step sizes to 
//' convolve. This typically incurs the least variance inflation. It is later 
//' also referred to as \emph{recursive nearest neighbor order of convolution}.
//' 
//' @inherit sequentialConvAll
//' 
//' @example inst/examples/closestPairConvAll.R
//'
// [[Rcpp::export]]
DataFrame closestPairConvAll(List lisOfDists, 
                             String convMethod = "splitAtom4productsIrregular", 
                             double headTrunc = 1e-10, 
                             double tailTrunc = 1e-10,
                             String regridMethod = "lr", 
                             unsigned N = 256,
                             unsigned forDeltaCompare = 256, 
                             bool useFFT = 0) 
{
  
  
  checkPMFlistIntegrity(lisOfDists);
  
  
  bool single = 1;
  
  if (lisOfDists.size() == 1) {
    List rst = lisOfDists[0];
    return rst;
  }

  std::vector<dist> lisOfDistsCpp(lisOfDists.size());
  int k = 0;
  double degenerateValSum = 0;
  for (int i = 0; i != lisOfDists.size(); ++i) {
    List tmp = lisOfDists[i];
    NumericVector value = tmp[0], prob = tmp[1];

    if (value.size() > 2) {
      lisOfDistsCpp[k].val.assign(value.begin(), value.end());
      lisOfDistsCpp[k].P.assign(prob.begin(), prob.end());
      k += 1;
    } 
    else if (value.size() == 2)
    {
      lisOfDistsCpp[k].val.resize(3);
      lisOfDistsCpp[k].val[0] = value[0];
      lisOfDistsCpp[k].val[1] = (value[0] + value[1]) / 2;
      lisOfDistsCpp[k].val[2] = value[1];
      lisOfDistsCpp[k].P.resize(3);
      lisOfDistsCpp[k].P[0] = prob[0];
      lisOfDistsCpp[k].P[1] = 0;
      lisOfDistsCpp[k].P[2] = prob[1];
      k += 1;
    }
    else
    {
      degenerateValSum += value[0];
    }
  }
  lisOfDistsCpp.resize(k);
  
  
  std::vector<pointerToDist> order(lisOfDistsCpp.size());
  std::vector<dist>::iterator i = lisOfDistsCpp.begin();
  std::vector<pointerToDist>::iterator j = order.begin();
  for (; i != lisOfDistsCpp.end(); ++i) {
    if (i->val.size() > 1) {
      j->range = i->val.back() - i->val.front();
      j->ptr = i;
      ++j;
    }
  }
  order.resize(j - order.begin());
  std::sort(order.begin(), order.end());
  dist result;
  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);
  closestPairConv(order, convMethod, regridMethod, N, forDeltaCompare, single,
                  headTrunc, tailTrunc, result, w);
  
  if (useFFT) delete w;

  
  double sumP = std::accumulate(result.P.begin(), result.P.end(), 0.0);
  for (std::vector<probtype>::iterator i = result.P.begin();
       i != result.P.end(); ++i)
    *i /= sumP;
  
  
  NumericVector finalRst(result.val.size());
  for (int i = 0, iend = result.val.size(); i < iend; ++i)
    finalRst[i] = result.val[i] + degenerateValSum;
  

  return DataFrame::create(Named("val") = finalRst, Named("P") = result.P);
}


double meanOfDist(dist &x) 
{
  unsigned end = x.val.size();
  double m = 0;
  for (unsigned i = 0; i != end; ++i) m += x.val[i] * x.P[i];
  return m;
}




// 3-way partition algorithm
dist partition3wayLRMconvCpp(std::vector<pointerToDist> &order,
                             String &convolutionMethod, String &regridMethod,
                             unsigned N, unsigned forDeltaCompare, bool single,
                             double headTrunc, double tailTrunc, wToUseFFT *w) 
{
  
  
  convMethod conv;
  
  
  if (convolutionMethod == "splitAtom4products")
    conv = &conv4productsCppTouchAtom;
  else if (convolutionMethod == "splitAtom9products")
    conv = &conv9productsCppTouchAtom;
  else if (convolutionMethod == "splitAtom4productsAtomUntouched")
    conv = &conv4productsCppAtomUntouched;
  else if (convolutionMethod == "splitAtom9productsAtomUntouched")
    conv = &conv9productsCppAtomUntouched;
  else if (convolutionMethod == "splitAtom4productsIrregular")
    conv = &conv4productsIrregularCpp;
  else if (convolutionMethod == "splitAtom9productsIrregular")
    conv = &conv9productsIrregularCpp;
  else conv = &convBtRgdForTest;

  
  truncMethod truncate;
  if (single) truncate = &truncateDistSingle;
  else truncate = &truncateDist;

  
  rgMethod rglinear_add = nullptr;
  if      (regridMethod == "lr")  rglinear_add = &rglinearAddCpp;
  else if (regridMethod == "r4")  rglinear_add = &rglr4splitAddCpp;
  else if (regridMethod == "lmm") rglinear_add = &lmmAddCpp;
  

  if (order.size() == 2) {
    dist tmp;
    tmp.val.resize(N);
    tmp.P.resize(N);
    tmp.getGrid(order[0].ptr->val.front() + order[1].ptr->val.front(),
                order[0].ptr->val.back() + order[1].ptr->val.back());

    conv(order[0].ptr->val.begin(), order[0].ptr->P.begin(),
         order[0].ptr->P.size(), order[1].ptr->val.begin(),
         order[1].ptr->P.begin(), order[1].ptr->P.size(), tmp.val.begin(),
         tmp.P.begin(), N, forDeltaCompare, tmp.val, tmp.P, rglinear_add, w);

    truncate(tmp.val, tmp.P, headTrunc, tailTrunc);
    return tmp;
  }

  
  // The first round
  std::vector<pointerToDist> middlePart(
      order.begin() + unsigned(order.size() / 3),
      order.begin() + unsigned(order.size() / 3) * 2);
  

  std::vector<pointerToDist>::iterator 
    leftBegin = order.begin(),
      leftEnd = order.begin() + unsigned(order.size() / 3),
      rightBegin = leftEnd + unsigned(order.size() / 3),
      rightEnd = order.end() - 1, compensate;

  
  unsigned sizeOfTuples = order.size();
  double sumOfMean, sumReserve = 0;

  
  for (std::vector<pointerToDist>::iterator i = order.begin(); 
       i != order.end(); ++i)
    sumReserve += i->range;
  sumOfMean = sumReserve / order.size() * 3;

  
  dist tmp(N);

  
  std::vector<dist> convStore1(
      order.size() - unsigned(order.size() / 3) * 2, dist(N)), 
      &convStore = convStore1;
  std::vector<dist>::iterator convStorei = convStore.begin();


  while (true) 
  {
    
    double forClosestToMean = sumOfMean - rightEnd->range - leftBegin->range;

    
    compensate = std::lower_bound(middlePart.begin(), middlePart.end(),
                                  forClosestToMean);

    if ((compensate != middlePart.begin() &&
         (compensate->range - forClosestToMean >
          forClosestToMean - (compensate - 1)->range)) ||
        compensate == middlePart.end())
      --compensate;
    

    tmp.val.resize(N);
    tmp.P.resize(N);
    tmp.getGrid(leftBegin->ptr->val.front() + compensate->ptr->val.front(),
                leftBegin->ptr->val.back() + compensate->ptr->val.back());

    conv(leftBegin->ptr->val.begin(), leftBegin->ptr->P.begin(),
         leftBegin->ptr->P.size(), compensate->ptr->val.begin(),
         compensate->ptr->P.begin(), compensate->ptr->P.size(), tmp.val.begin(),
         tmp.P.begin(), N, forDeltaCompare, tmp.val, tmp.P, rglinear_add, w);

    truncate(tmp.val, tmp.P, headTrunc, tailTrunc);


    convStorei->val.resize(N);
    convStorei->P.resize(N);
    convStorei->getGrid(tmp.val.front() + rightEnd->ptr->val.front(),
                        tmp.val.back() + rightEnd->ptr->val.back());

    conv(tmp.val.begin(), tmp.P.begin(), tmp.P.size(),
         rightEnd->ptr->val.begin(), rightEnd->ptr->P.begin(),
         rightEnd->ptr->P.size(), convStorei->val.begin(),
         convStorei->P.begin(), N, forDeltaCompare, convStorei->val,
         convStorei->P, rglinear_add, w);

    truncate(convStorei->val, convStorei->P, headTrunc, tailTrunc);
    
    
    // replce the range of order
    leftBegin->range = convStorei->val.back() - convStorei->val.front();
    leftBegin->ptr = convStorei;


    ++convStorei;
    ++leftBegin;
    --rightEnd;
    if (leftBegin == leftEnd) break;
    // left part is always smaller than or equal to the right part
    middlePart.erase(compensate);
  }
  

  for (; rightEnd >= rightBegin; ++convStorei, --rightEnd, ++leftBegin) {
    convStorei->val.assign(rightEnd->ptr->val.begin(),
                           rightEnd->ptr->val.end());
    convStorei->P.assign(rightEnd->ptr->P.begin(), rightEnd->ptr->P.end());
    leftBegin->ptr = convStorei;
    leftBegin->range = convStorei->val.back() - convStorei->val.front();
  }

  
  std::vector<dist> convStore2(
      convStore.size() - unsigned(convStore.size() / 3) * 2, dist(N)),
      &convStoreAlter = convStore2;
  std::vector<dist>::iterator convStoreAlteri;

  
  while (true) 
  {
    sizeOfTuples = convStore.size();
  
    if (sizeOfTuples == 1 || sizeOfTuples == 2) break;

    order.resize(sizeOfTuples);
    convStoreAlter.resize(sizeOfTuples - unsigned(sizeOfTuples / 3) * 2);

    sort(order.begin(), order.end());
    // middlePart won't be reallocated
    middlePart.assign(order.begin() + unsigned(convStore.size() / 3),
                      order.begin() + unsigned(convStore.size() / 3) * 2);

    sumOfMean = sumReserve / sizeOfTuples * 3;

    leftBegin = order.begin();
    leftEnd = leftBegin + unsigned(convStore.size() / 3);
    
    
    rightBegin = leftEnd + unsigned(convStore.size() / 3);
    rightEnd = order.end() - 1;
    convStoreAlteri = convStoreAlter.begin();


    while (true) {
      double forClosestToMean = sumOfMean - rightEnd->range - leftBegin->range;
      compensate = std::lower_bound(middlePart.begin(), middlePart.end(),
                                    forClosestToMean);

      if ((compensate != middlePart.begin() &&
           (compensate->range - forClosestToMean >
            forClosestToMean - (compensate - 1)->range)) ||
          compensate == middlePart.end())
        --compensate;
      // convolve the distributions pointer by leftBegin, compensate and
      // rightEnd..

      // myfile<<"compensate
      // position=="<<compensate-middlePart.begin()<<std::endl;
      // myfile<<"leftBegin position=="<<leftBegin-order.begin()<<std::endl;
      // myfile<<"rightEnd position=="<<rightEnd-(order.end()-1)<<std::endl;

      tmp.val.resize(N);
      tmp.P.resize(N);

      tmp.getGrid(leftBegin->ptr->val.front() + compensate->ptr->val.front(),
                  leftBegin->ptr->val.back() + compensate->ptr->val.back());

      conv(leftBegin->ptr->val.begin(), leftBegin->ptr->P.begin(),
           leftBegin->ptr->P.size(), compensate->ptr->val.begin(),
           compensate->ptr->P.begin(), compensate->ptr->P.size(),
           tmp.val.begin(), tmp.P.begin(), N, forDeltaCompare, tmp.val, tmp.P,
           rglinear_add, w);

      truncate(tmp.val, tmp.P, headTrunc, tailTrunc);

      convStoreAlteri->val.resize(N);
      convStoreAlteri->P.resize(N);

      convStoreAlteri->getGrid(tmp.val.front() + rightEnd->ptr->val.front(),
                               tmp.val.back() + rightEnd->ptr->val.back());

      conv(tmp.val.begin(), tmp.P.begin(), tmp.P.size(),
           rightEnd->ptr->val.begin(), rightEnd->ptr->P.begin(),
           rightEnd->ptr->P.size(), convStoreAlteri->val.begin(),
           convStoreAlteri->P.begin(), N, forDeltaCompare, convStoreAlteri->val,
           convStoreAlteri->P, rglinear_add, w);

      truncate(convStoreAlteri->val, convStoreAlteri->P, headTrunc, tailTrunc);

      // leftBegin->range+=rightEnd->range+compensate->range;
      leftBegin->ptr = convStoreAlteri;
      leftBegin->range =
          convStoreAlteri->val.back() - convStoreAlteri->val.front();

      ++leftBegin;
      --rightEnd;
      ++convStoreAlteri;

      if (leftBegin == leftEnd) break;
      // left part is always smaller than or equal to the right part
      middlePart.erase(compensate);
    }


    for (; rightEnd >= rightBegin; ++convStoreAlteri, --rightEnd, ++leftBegin) 
    {
      convStoreAlteri->val.assign(rightEnd->ptr->val.begin(),
                                  rightEnd->ptr->val.end());
      convStoreAlteri->P.assign(rightEnd->ptr->P.begin(),
                                rightEnd->ptr->P.end());
      leftBegin->ptr = convStoreAlteri;
      leftBegin->range =
          convStoreAlteri->val.back() - convStoreAlteri->val.front();
    }

    convStore.resize(convStoreAlter.size());
    std::swap(convStore, convStoreAlter);
  }


  if (sizeOfTuples == 2) 
  {
    tmp.val.resize(N);
    tmp.P.resize(N);
    tmp.getGrid(convStore.front().val.front() + convStore.back().val.front(),
                convStore.front().val.back() + convStore.back().val.back());

    
    std::fill(tmp.P.begin(), tmp.P.end(), 0);

    
    conv(convStore.front().val.begin(), convStore.front().P.begin(),
         convStore.front().P.size(), convStore.back().val.begin(),
         convStore.back().P.begin(), convStore.back().P.size(), tmp.val.begin(),
         tmp.P.begin(), N, forDeltaCompare, tmp.val, tmp.P, rglinear_add, w);

    truncate(tmp.val, tmp.P, headTrunc, tailTrunc);
    return tmp;
  }
  return convStore.front();
}




// [[Rcpp::export]]
DataFrame partition3wayLRMconv(List lisOfDists, String convMethod, bool single,
                               double headTrunc, double tailTrunc,
                               String regridMethod = "lr", unsigned N = 256,
                               unsigned forDeltaCompare = 256,
                               bool useFFT = 0) 
{
  if (lisOfDists.size() == 1) {
    List rst = lisOfDists[0];
    return rst;
  }

  std::vector<dist> lisOfDistsCpp(lisOfDists.size());
  for (int i = 0; i != lisOfDists.size(); ++i) {
    List tmp = lisOfDists[i];
    NumericVector value = tmp[0], prob = tmp[1];

    if (value.size() > 2) {
      lisOfDistsCpp[i].val.assign(value.begin(), value.end());
      lisOfDistsCpp[i].P.assign(prob.begin(), prob.end());
    } else {
      lisOfDistsCpp[i].val.resize(3);
      lisOfDistsCpp[i].val[0] = value[0];
      lisOfDistsCpp[i].val[1] = (value[0] + value[1]) / 2;
      lisOfDistsCpp[i].val[2] = value[1];
      lisOfDistsCpp[i].P.resize(3);
      lisOfDistsCpp[i].P[0] = prob[0];
      lisOfDistsCpp[i].P[1] = 0;
      lisOfDistsCpp[i].P[2] = prob[1];
    }
  }
  // std::vector<double>rst(N), rstp(N,0);

  std::vector<pointerToDist> order(lisOfDistsCpp.size());
  std::vector<dist>::iterator i = lisOfDistsCpp.begin();
  std::vector<pointerToDist>::iterator j = order.begin();
  for (; i != lisOfDistsCpp.end(); ++i) {
    if (i->val.size() > 1) {
      j->range = i->val.back() - i->val.front();
      j->ptr = i;
      ++j;
    }
  }
  order.resize(j - order.begin());
  std::sort(order.begin(), order.end());

  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);

  dist result =
      partition3wayLRMconvCpp(order, convMethod, regridMethod, N,
                              forDeltaCompare, single, headTrunc, tailTrunc, w);
  delete w;

  double sumP = std::accumulate(result.P.begin(), result.P.end(), 0.0);
  for (std::vector<probtype>::iterator i = result.P.begin();
       i != result.P.end(); ++i)
    *i /= sumP;
  return DataFrame::create(Named("val") = result.val, Named("P") = result.P);
}




dist BLDMconvCpp(std::vector<pointerToDist> &order, /*std::vector<double>&rst,
                     std::vector<double>&rstp,*/
                 String &convolutionMethod, String &regridMethod, unsigned N,
                 unsigned forDeltaCompare, bool single, double headTrunc,
                 double tailTrunc, wToUseFFT *w) 
{
  convMethod conv;
  if (convolutionMethod == "splitAtom4products")
    conv = &conv4productsCppTouchAtom;
  else if (convolutionMethod == "splitAtom9products")
    conv = &conv9productsCppTouchAtom;
  else if (convolutionMethod == "splitAtom4productsAtomUntouched")
    conv = &conv4productsCppAtomUntouched;
  else if (convolutionMethod == "splitAtom9productsAtomUntouched")
    conv = &conv9productsCppAtomUntouched;
  else if (convolutionMethod == "splitAtom4productsIrregular")
    conv = &conv4productsIrregularCpp;
  else if (convolutionMethod == "splitAtom9productsIrregular")
    conv = &conv9productsIrregularCpp;
  else
    conv = &convBtRgdForTest;

  truncMethod truncate;
  if (single)
    truncate = &truncateDistSingle;
  else
    truncate = &truncateDist;

  rgMethod rglinear_add = nullptr;
  if (regridMethod == "lr") rglinear_add = &rglinearAddCpp;
  // else if(regridMethod=="r3")rglinear_add=&rglr3splitAddCpp;
  else if (regridMethod == "r4")
    rglinear_add = &rglr4splitAddCpp;
  else if (regridMethod == "lmm")
    rglinear_add = &lmmAddCpp;
  // else rglinear_add=&rglr4splitGaussianKernalAddCpp;

  if (order.size() == 2) {
    dist tmp;
    tmp.val.resize(N);
    tmp.P.resize(N);
    tmp.getGrid(order[0].ptr->val.front() + order[1].ptr->val.front(),
                order[0].ptr->val.back() + order[1].ptr->val.back());

    conv(order[0].ptr->val.begin(), order[0].ptr->P.begin(),
         order[0].ptr->P.size(), order[1].ptr->val.begin(),
         order[1].ptr->P.begin(), order[1].ptr->P.size(), tmp.val.begin(),
         tmp.P.begin(), N, forDeltaCompare, tmp.val, tmp.P, rglinear_add, w);

    truncate(tmp.val, tmp.P, headTrunc, tailTrunc);
    return tmp;
  }

  std::vector<dist> convStore1(order.size() - unsigned(order.size() / 2),
                               dist(N)),
      &convStore = convStore1,
      convStore2(convStore1.size() - unsigned(convStore1.size() / 2), dist(N)),
      &convStoreAlter = convStore2;

  std::vector<pointerToDist>::iterator left = order.begin(),
                                       right = order.end() - 1;
  std::vector<dist>::iterator convStorei = convStore.begin();
  for (; left < right; ++left, --right, ++convStorei) {
    convStorei->val.resize(N);
    convStorei->P.resize(N);
    convStorei->getGrid(left->ptr->val.front() + right->ptr->val.front(),
                        left->ptr->val.back() + right->ptr->val.back());
    conv(left->ptr->val.begin(), left->ptr->P.begin(), left->ptr->P.size(),
         right->ptr->val.begin(), right->ptr->P.begin(), right->ptr->P.size(),
         convStorei->val.begin(), convStorei->P.begin(), N, forDeltaCompare,
         convStorei->val, convStorei->P, rglinear_add, w);
    truncate(convStorei->val, convStorei->P, headTrunc, tailTrunc);
    left->ptr = convStorei;
    left->range = convStorei->val.back() - convStorei->val.front();
  }
  if (left == right) {
    convStorei->val.assign(left->ptr->val.begin(), left->ptr->val.end());
    convStorei->P.assign(left->ptr->P.begin(), left->ptr->P.end());
    left->ptr = convStorei;
    left->range = convStorei->val.back() - convStorei->val.front();
  }

  std::vector<dist>::iterator convStoreAlteri;

  while (true) 
  {
    order.resize(convStore.size());
    std::sort(order.begin(), order.end());
    left = order.begin();
    right = order.end() - 1;
    convStoreAlteri = convStoreAlter.begin();
    for (; left < right; ++left, --right, ++convStoreAlteri) 
    {
      convStoreAlteri->val.resize(N);
      convStoreAlteri->P.resize(N);
      convStoreAlteri->getGrid(left->ptr->val.front() + right->ptr->val.front(),
                               left->ptr->val.back() + right->ptr->val.back());
      conv(left->ptr->val.begin(), left->ptr->P.begin(), left->ptr->P.size(),
           right->ptr->val.begin(), right->ptr->P.begin(), right->ptr->P.size(),
           convStoreAlteri->val.begin(), convStoreAlteri->P.begin(), N,
           forDeltaCompare, convStoreAlteri->val, convStoreAlteri->P,
           rglinear_add, w);
      truncate(convStoreAlteri->val, convStoreAlteri->P, headTrunc, tailTrunc);
      left->ptr = convStoreAlteri;
      left->range = convStoreAlteri->val.back() - convStoreAlteri->val.front();
    }
    if (convStoreAlter.size() == 1) break;
    if (left == right) 
    {
      convStoreAlteri->val.assign(left->ptr->val.begin(), left->ptr->val.end());
      convStoreAlteri->P.assign(left->ptr->P.begin(), left->ptr->P.end());
      left->ptr = convStoreAlteri;
      left->range = convStoreAlteri->val.back() - convStoreAlteri->val.front();
    }
    std::swap(convStore, convStoreAlter);
    convStoreAlter.resize(convStore.size() - unsigned(convStore.size() / 2));
  }

  
  return convStoreAlter.front();
}


// [[Rcpp::export]]
DataFrame BLDMconv(List lisOfDists, String convMethod, bool single,
                   double headTrunc, double tailTrunc,
                   String regridMethod = "lr", unsigned N = 256,
                   unsigned forDeltaCompare = 256, bool useFFT = 0) {
  if (lisOfDists.size() == 1) {
    List rst = lisOfDists[0];
    return rst;
  }

  std::vector<dist> lisOfDistsCpp(lisOfDists.size());
  for (int i = 0; i != lisOfDists.size(); ++i) {
    List tmp = lisOfDists[i];
    NumericVector value = tmp[0], prob = tmp[1];

    if (value.size() > 2) {
      lisOfDistsCpp[i].val.assign(value.begin(), value.end());
      lisOfDistsCpp[i].P.assign(prob.begin(), prob.end());
    } else {
      lisOfDistsCpp[i].val.resize(3);
      lisOfDistsCpp[i].val[0] = value[0];
      lisOfDistsCpp[i].val[1] = (value[0] + value[1]) / 2;
      lisOfDistsCpp[i].val[2] = value[1];
      lisOfDistsCpp[i].P.resize(3);
      lisOfDistsCpp[i].P[0] = prob[0];
      lisOfDistsCpp[i].P[1] = 0;
      lisOfDistsCpp[i].P[2] = prob[1];
    }
  }

  std::vector<pointerToDist> order(lisOfDistsCpp.size());
  std::vector<dist>::iterator i = lisOfDistsCpp.begin();
  std::vector<pointerToDist>::iterator j = order.begin();
  for (; i != lisOfDistsCpp.end(); ++i) {
    if (i->val.size() > 1) {
      j->range = i->val.back() - i->val.front();
      j->ptr = i;
      ++j;
    }
  }
  order.resize(j - order.begin());
  std::sort(order.begin(), order.end());

  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);

  dist result = BLDMconvCpp(order, convMethod, regridMethod, N, forDeltaCompare,
                            single, headTrunc, tailTrunc, w);
  delete w;
  double sumP = std::accumulate(result.P.begin(), result.P.end(), 0.0);
  for (std::vector<probtype>::iterator i = result.P.begin();
       i != result.P.end(); ++i)
    *i /= sumP;

  // Rcout<<lrCount<<"\n";

  return DataFrame::create(Named("val") = result.val, Named("P") = result.P);
}

dist forBiggestRangeDiffCpp(std::vector<pointerToDist> &order,
                            /*std::vector<double>&rst,
std::vector<double>&rstp,*/ String &convolutionMethod, String &regridMethod,
                            unsigned N, unsigned forDeltaCompare, bool single,
                            double headTrunc, double tailTrunc, wToUseFFT *w) {
  convMethod conv;
  if (convolutionMethod == "splitAtom4products")
    conv = &conv4productsCppTouchAtom;
  else if (convolutionMethod == "splitAtom9products")
    conv = &conv9productsCppTouchAtom;
  else if (convolutionMethod == "splitAtom4productsAtomUntouched")
    conv = &conv4productsCppAtomUntouched;
  else if (convolutionMethod == "splitAtom9productsAtomUntouched")
    conv = &conv9productsCppAtomUntouched;
  else if (convolutionMethod == "splitAtom4productsIrregular")
    conv = &conv4productsIrregularCpp;
  else if (convolutionMethod == "splitAtom9productsIrregular")
    conv = &conv9productsIrregularCpp;
  else
    conv = &convBtRgdForTest;

  truncMethod truncate;
  if (single)
    truncate = &truncateDistSingle;
  else
    truncate = &truncateDist;

  rgMethod rglinear_add = nullptr;
  if (regridMethod == "lr") rglinear_add = &rglinearAddCpp;
  // else if(regridMethod=="r3")rglinear_add=&rglr3splitAddCpp;
  else if (regridMethod == "r4")
    rglinear_add = &rglr4splitAddCpp;
  else if (regridMethod == "lmm")
    rglinear_add = &lmmAddCpp;
  // else rglinear_add=&rglr4splitGaussianKernalAddCpp;

  if (order.size() == 2) {
    dist tmp;
    tmp.val.resize(N);
    tmp.P.resize(N);
    tmp.getGrid(order[0].ptr->val.front() + order[1].ptr->val.front(),
                order[0].ptr->val.back() + order[1].ptr->val.back());

    conv(order[0].ptr->val.begin(), order[0].ptr->P.begin(),
         order[0].ptr->P.size(), order[1].ptr->val.begin(),
         order[1].ptr->P.begin(), order[1].ptr->P.size(), tmp.val.begin(),
         tmp.P.begin(), N, forDeltaCompare, tmp.val, tmp.P, rglinear_add, w);

    truncate(tmp.val, tmp.P, headTrunc, tailTrunc);
    return tmp;
  }

  dist tmp(N), &convStore = tmp, tmp2(N), &convStoreAlter = tmp2;

  tmp.val.resize(N);
  tmp.P.resize(N);
  tmp.getGrid(order.front().ptr->val.front() + order.back().ptr->val.front(),
              order.front().ptr->val.back() + order.back().ptr->val.back());

  conv(order.front().ptr->val.begin(), order.front().ptr->P.begin(),
       order.front().ptr->P.size(), order.back().ptr->val.begin(),
       order.back().ptr->P.begin(), order.back().ptr->P.size(), tmp.val.begin(),
       tmp.P.begin(), N, forDeltaCompare, tmp.val, tmp.P, rglinear_add, w);
  truncate(tmp.val, tmp.P, headTrunc, tailTrunc);

  std::vector<pointerToDist>::iterator i = order.begin() + 1,
                                       end = order.end() - 1;

  for (; i < end; ++i) {
    convStoreAlter.val.resize(N);
    convStoreAlter.P.resize(N);
    convStoreAlter.getGrid(convStore.val.front() + i->ptr->val.front(),
                           convStore.val.back() + i->ptr->val.back());
    conv(convStore.val.begin(), convStore.P.begin(), convStore.P.size(),
         i->ptr->val.begin(), i->ptr->P.begin(), i->ptr->P.size(),
         convStoreAlter.val.begin(), convStoreAlter.P.begin(), N,
         forDeltaCompare, convStoreAlter.val, convStoreAlter.P, rglinear_add,
         w);
    truncate(tmp.val, tmp.P, headTrunc, tailTrunc);
    std::swap(convStore, convStoreAlter);
  }
  return convStore;
}

// [[Rcpp::export]]
DataFrame forBiggestRangeDiff(List lisOfDists, String convMethod, bool single,
                              double headTrunc, double tailTrunc,
                              String regridMethod = "lr", unsigned N = 256,
                              unsigned forDeltaCompare = 256, bool useFFT = 0) {
  if (lisOfDists.size() == 1) {
    List rst = lisOfDists[0];
    return rst;
  }

  std::vector<dist> lisOfDistsCpp(lisOfDists.size());
  for (int i = 0; i != lisOfDists.size(); ++i) {
    List tmp = lisOfDists[i];
    NumericVector value = tmp[0], prob = tmp[1];

    if (value.size() > 2) {
      lisOfDistsCpp[i].val.assign(value.begin(), value.end());
      lisOfDistsCpp[i].P.assign(prob.begin(), prob.end());
    } else {
      lisOfDistsCpp[i].val.resize(3);
      lisOfDistsCpp[i].val[0] = value[0];
      lisOfDistsCpp[i].val[1] = (value[0] + value[1]) / 2;
      lisOfDistsCpp[i].val[2] = value[1];
      lisOfDistsCpp[i].P.resize(3);
      lisOfDistsCpp[i].P[0] = prob[0];
      lisOfDistsCpp[i].P[1] = 0;
      lisOfDistsCpp[i].P[2] = prob[1];
    }
  }

  std::vector<pointerToDist> order(lisOfDistsCpp.size());
  std::vector<dist>::iterator i = lisOfDistsCpp.begin();
  std::vector<pointerToDist>::iterator j = order.begin();
  for (; i != lisOfDistsCpp.end(); ++i) 
  {
    if (i->val.size() > 1) {
      j->range = i->val.back() - i->val.front();
      j->ptr = i;
      ++j;
    }
  }
  order.resize(j - order.begin());
  std::sort(order.begin(), order.end());
  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);
  dist result =
      forBiggestRangeDiffCpp(order, convMethod, regridMethod, N,
                             forDeltaCompare, single, headTrunc, tailTrunc, w);
  delete w;
  double sumP = std::accumulate(result.P.begin(), result.P.end(), 0.0);
  for (std::vector<probtype>::iterator i = result.P.begin();
       i != result.P.end(); ++i)
    *i /= sumP;
  return DataFrame::create(Named("val") = result.val, Named("P") = result.P);
}


//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
//---------------------------mix-convolution with spatial correlation:


// [[Rcpp::export]]
IntegerVector lowerBound(NumericVector x, NumericVector v) 
{
  IntegerVector rst(x.size());
  for (int i = 0, iend = x.size(); i < iend; ++i) {
    rst[i] = std::lower_bound(v.begin(), v.end(), x[i]) - v.begin() + 1;
  }
  return rst;
}


double Mean(std::vector<losstype> &val, std::vector<probtype> &P) 
{
  return std::inner_product(val.begin(), val.end(), P.begin(), 0.0);
}


double Var(std::vector<losstype> &val, std::vector<probtype> &P) {
  double S = 0, m = Mean(val, P);
  for (unsigned i = 0, iend = val.size(); i != iend; ++i)
    S += (val[i] * val[i] * P[i]);
  return S - m * m;
}



void comoPtr(losstype *x, probtype *xp, int xsize,
             losstype *y, probtype *yp, int ysize,
             std::vector<losstype> &val, std::vector<probtype> &P)
{
  // rst and rstp should at first have the size of at least xlen+ylen
  int xi = 0, yi = 0, rei = 0;
  // unsigned xsize = x.size(), ysize = y.size();
  losstype tmp;
  while (true) {
    val[rei] = x[xi] + y[yi];
    tmp = xp[xi] - yp[yi];
    
    if (tmp > 0) {
      P[rei] = yp[yi];
      xp[xi] = tmp;
      ++yi;
      
      if (yi == ysize) {
        {  // tackle the numeric error in case there exists any
          if (xi < xsize - 1)  // meanwhile x pointer haven't reached the tail
          {
            for (xi += 1, rei += 1;;) {
              val[rei] = x[xi] + y[ysize - 1];
              P[rei] = xp[xi];
              ++xi;
              if (xi >= xsize)
                break;  // don't change the structure. We are making sure rei
              // point to the last element, not the end
              ++rei;
            }
          }
        }
        break;
      }
      
    } else if (tmp < 0) {
      P[rei] = xp[xi];
      yp[yi] = -tmp;
      ++xi;
      
      if (xi == xsize) {
        {  // tackle the numeric error in case there exists any
          if (yi < ysize - 1)  // meanwhile y pointer haven't reached the tail
          {
            for (yi += 1, ++rei;;) {
              val[rei] = y[yi] + x[xsize - 1];
              P[rei] = yp[yi];
              ++yi;
              if (yi >= ysize) break;
              ++rei;
            }
          }
        }
        
        break;
      }
      
    } else {
      P[rei] = xp[xi];
      ++xi;
      ++yi;
      
      if (xi == xsize) {
        if (yi < ysize - 1)  // meanwhile y pointer haven't reached the tail
        {
          for (yi += 1, ++rei;;) {
            val[rei] = y[yi] + x[xsize - 1];
            P[rei] = yp[yi];
            ++yi;
            if (yi >= ysize) break;
            ++rei;
          }
        }
        
        break;
      }
      
      if (yi == ysize) {
        if (xi < xsize - 1)  // meanwhile x pointer haven't reached the tail
        {
          for (xi += 1, ++rei;;) {
            val[rei] = x[xi] + y[ysize - 1];
            P[rei] = xp[xi];
            ++xi;
            if (xi >= xsize) break;
            ++rei;
          }
        }
        
        break;
      }
    }
    ++rei;
  }
  
  // sometimes the sum of probabilities are strictly equal to 1 for both
  // distribution because numeric error exists. This can result in the
  // comonotonization's max loss point does not equal to the sum of the
  // individual distributions' max loss points
  
  val.resize(rei + 1);
  P.resize(rei + 1);
}


void comoVec(std::vector<losstype> &x, std::vector<probtype> &xp,
          std::vector<losstype> &y, std::vector<probtype> &yp,
          std::vector<losstype> &val, std::vector<probtype> &P) 
{
  comoPtr(&x[0], &xp[0], x.size(), &y[0], &yp[0], y.size(), val, P);
}




//'
//' Comonotonization
//' 
//' Compute the sum of two distributions assuming maximum correlation.
//' 
//' @inheritParams convSplitAtom4productsIrregular
//' 
//' @param N  Support size of the result PMF. Default 64. If \code{<= 1}, 
//' the function comonotonizes \code{X} and \code{Y} and returns the 
//' result without regriding. Comonotonization does not need the two PMFs to 
//' have common step sizes. The regriding can prevent the result PMF from being
//' too long.
//'
//' @param rgMethod  Regriding method:
//' 
//'   \code{"lr"}: linear regrid. Default.
//' 
//'   \code{"r4"}: four-point regrid.
//' 
//'   \code{""}: do not regrid. Has the same effect as \code{N <= 1}.
//' 
//' @param isUsedInCounterComo TRUE if this procedure is used for
//' counter-comonotonization. Default FALSE.
//'
//' @inherit rglr return
//' 
//' @inherit convSplitAtom4productsIrregular details
//' 
//' @example inst/examples/como.R
//' 
// [[Rcpp::export]]
DataFrame como(DataFrame X, DataFrame Y, int N = 64, String rgMethod = "lr",
               bool isUsedInCounterComo = false) 
{
  
  
  checkPMFintegrity(X, isUsedInCounterComo);
  checkPMFintegrity(Y, isUsedInCounterComo);
  
  
  NumericVector x = X[0], y = Y[0];
  NumericVector xp_ = X[1], yp_ = Y[1];
  std::vector<probtype> xp(xp_.begin(), xp_.end()), yp(yp_.begin(), yp_.end());
  
  
  std::vector<losstype> rst(x.size() + y.size());
  std::vector<probtype> rstp(x.size() + y.size());
  comoPtr(&x[0], &xp[0], x.size(), &y[0], &yp[0], y.size(), rst, rstp);
  if (N <= 1 or rgMethod == "") 
    return DataFrame::create(Named("val") = rst, Named("P") = rstp);
  
  
  NumericVector ngdP(N), ngd(N);
  double ngdBegin = x[0] + y[0];
  double span = *(x.end() - 1) + *(y.end() - 1) - ngdBegin, 
    delta = span / (N - 1);
  for (int i = 0; i < N; ++i) ngd[i] = ngdBegin + i * delta;
  
  
  if (rgMethod == "lr") 
  {
    rglinearAddCppPtr(&*rst.begin(), &*rstp.begin(), rst.size(), 
                      &*ngd.begin(), &*ngdP.begin(), ngd.size());
  } 
  else 
  {
    rg4(&rst[0], &rstp[0], rst.size(), &ngd[0], &ngdP[0], ngd.size());
  }
  return DataFrame::create(Named("val") = ngd, Named("P") = ngdP);
}


void comoJoint(std::vector<losstype> &x, std::vector<probtype> &xp,
               std::vector<losstype> &y, std::vector<probtype> &yp,
               // std::vector<losstype> &val,
               std::vector<std::pair<losstype, losstype> > &val,
               std::vector<probtype> &P) {
  // rst and rstp should at first have the size of at least xlen+ylen
  unsigned xi = 0, yi = 0, rei = 0, xsize = x.size(), ysize = y.size();
  losstype tmp;
  while (true) {
    // val[rei]=x[xi]+y[yi];
    val[rei] = std::pair<losstype, losstype>(x[xi], y[yi]);
    tmp = xp[xi] - yp[yi];

    if (tmp > 0) {
      P[rei] = yp[yi];
      xp[xi] = tmp;
      ++yi;

      if (yi == ysize) {
        {  // tackle the numeric error in case there exists any
          if (xi < xsize - 1)  // meanwhile x pointer haven't reached the tail
          {
            for (xi += 1, rei += 1;;) {
              // val[rei] = x[xi] + y[ysize-1];
              val[rei] = std::pair<losstype, losstype>(x[xi], y[ysize - 1]);
              P[rei] = xp[xi];
              ++xi;
              if (xi >= xsize)
                break;  // don't change the structure. We are making sure rei
                        // point to the last element, not the end
              ++rei;
            }
          }
        }
        break;
      }

    } else if (tmp < 0) {
      P[rei] = xp[xi];
      yp[yi] = -tmp;
      ++xi;
      if (xi == xsize) {
        {  // tackle the numeric error in case there exists any
          if (yi < ysize - 1)  // meanwhile y pointer haven't reached the tail
          {
            for (yi += 1, ++rei;;) {
              // val[rei]=y[yi]+x[xsize-1];
              val[rei] = std::pair<losstype, losstype>(x[xsize - 1], y[yi]);
              P[rei] = yp[yi];
              ++yi;
              if (yi >= ysize) break;
              ++rei;
            }
          }
        }
        break;
      }
    } else {
      P[rei] = xp[xi];
      ++xi;
      ++yi;
      if (xi == xsize) {
        if (yi < ysize - 1)  // meanwhile y pointer haven't reached the tail
        {
          for (yi += 1, ++rei;;) {
            // val[rei] = y[yi] + x[xsize - 1];
            val[rei] = std::pair<losstype, losstype>(x[xsize - 1], y[yi]);
            P[rei] = yp[yi];
            ++yi;
            if (yi >= ysize) break;
            ++rei;
          }
        }

        break;
      }

      if (yi == ysize) {
        if (xi < xsize - 1)  // meanwhile x pointer haven't reached the tail
        {
          for (xi += 1, ++rei;;) {
            // val[rei]=x[xi]+y[ysize-1];
            val[rei] = std::pair<losstype, losstype>(x[xi], y[ysize - 1]);
            P[rei] = xp[xi];
            ++xi;
            if (xi >= xsize) break;
            ++rei;
          }
        }

        break;
      }
    }

    ++rei;
  }

  // sometimes the sum of probabilities are strictly equal to 1 for both
  // distribution because numeric error exists. This can result in the
  // comonotonization's max loss point does not equal to the sum of the
  // individual distributions' max loss points

  val.resize(rei + 1);
  P.resize(rei + 1);

  // val.resize(rei);
  // P.resize(rei);
}


// [[Rcpp::export]]
DataFrame comoJoint2D(List X, List Y) {
  std::vector<losstype> x = X[0], y = Y[0];
  std::vector<probtype> xp = X[1], yp = Y[1];
  std::vector<std::pair<double, double> > rst(x.size() + y.size());
  std::vector<probtype> rstp(x.size() + y.size());
  comoJoint(x, xp, y, yp, rst, rstp);
  NumericVector jx(rst.size()), jy(rst.size());
  for (int k = 0, kend = rst.size(); k < kend; ++k) {
    jx[k] = rst[k].first;
    jy[k] = rst[k].second;
  }
  return DataFrame::create(Named("x") = jx, Named("y") = jy, Named("P") = rstp);
}

// this function will change x and y
// return the mix-convolution's variance
double mixConv(std::vector<losstype> &x, std::vector<probtype> &xp, double varx,
               std::vector<losstype> &y, std::vector<probtype> &yp, double vary,
               std::vector<losstype> &rst, std::vector<probtype> &rstp,
               unsigned forDeltaCompare, rgMethod rglinear_add,
               double supposedVar, unsigned &outOfBound, bool pureConv,
               wToUseFFT *w) 
{
  
  
  if (x.size() == 1) {
    std::swap(rst, y);
    std::swap(rstp, yp);
    return 0;
  }
  if (y.size() == 1) {
    std::swap(rst, x);
    std::swap(rstp, xp);
    return 0;
  }

  bool varianceEnlarged = conv4productsIrregularCpp(
      x.begin(), xp.begin(), x.size(), y.begin(), yp.begin(), y.size(),
      rst.begin(), rstp.begin(), rst.size(), forDeltaCompare, rst, rstp,
      rglinear_add, w);
  if (pureConv) {
    if (!varianceEnlarged) return varx + vary;
    return Var(rst, rstp);
  }

  double varConv;
  if (varianceEnlarged)
    varConv = Var(rst, rstp);
  else
    varConv = varx + vary;

  // bool varianceImposed=0;

  if (supposedVar / varConv > 1.0001) {
    std::vector<probtype> Pnew(rstp.size());
    {
      std::vector<losstype> val(x.size() + y.size());
      std::vector<probtype> P(val.size());
      comoVec(x, xp, y, yp, val, P);
      rglinear_add(val.begin(), P.begin(), P.size(), rst.begin(), Pnew.begin(),
                   rst.size());
    }
    // now como is stored in rst and Pnew
    double varComo = Var(rst, Pnew),
           w = (supposedVar - varConv) / (varComo - varConv);

    if (w <= 0) {
      ++outOfBound;
      return varConv;
    } else if (w >= 1) {
      std::swap(rstp, Pnew);
      ++outOfBound;
      return varComo;
    } else {
      // varianceImposed=1;
      for (unsigned i = 0, iend = rst.size(); i != iend; ++i)
        rstp[i] = (1 - w) * rstp[i] + w * Pnew[i];
      return supposedVar;
    }
  } else
    return varConv;
}


// [[Rcpp::export]]
DataFrame mixConv2dist(List dist1, List dist2, double corWanted,
                       char regridMethod = 'r', int N = 256,
                       double dist1var = -1, double dist2var = -1,
                       bool useFFT = 0) {
  std::vector<losstype> x = dist1[0], y = dist2[0];
  std::vector<probtype> xp = dist1[1], yp = dist2[1];

  if (dist1var < 0) dist1var = Var(x, xp);
  if (dist2var < 0) dist2var = Var(y, yp);
  double supposedVar =
      dist1var + dist2var + 2 * corWanted * std::sqrt(dist1var * dist2var);
  rgMethod regM;
  if (regridMethod == 'l')
    regM = rglinearAddCpp;
  else
    regM = rglr4splitAddCpp;
  unsigned outOfBound = 0;
  std::vector<losstype> rst(N);
  std::vector<probtype> rstp(N);
  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);
  mixConv(x, xp, dist1var, y, yp, dist2var, rst, rstp, N, regM, supposedVar,
          outOfBound, 0, w);
  delete w;
  if (outOfBound != 0) std::cout << "weight out of bound\n";
  return DataFrame::create(Named("val") = rst, Named("P") = rstp);
}

// [[Rcpp::export]]
DataFrame mixConv2distGivenVar(List dist1, List dist2, double varWanted,
                               char regridMethod = 'r', int N = 256,
                               double dist1var = -1, double dist2var = -1,
                               bool useFFT = 0) {
  std::vector<losstype> x = dist1[0], y = dist2[0];
  std::vector<probtype> xp = dist1[1], yp = dist2[1];
  if (dist1var < 0) dist1var = Var(x, xp);
  if (dist2var < 0) dist2var = Var(y, yp);
  double supposedVar = varWanted;
  rgMethod regM;
  if (regridMethod == 'l')
    regM = rglinearAddCpp;
  else
    regM = rglr4splitAddCpp;
  unsigned outOfBound = 0;
  std::vector<losstype> rst(N);
  std::vector<probtype> rstp(N);
  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);
  mixConv(x, xp, dist1var, y, yp, dist2var, rst, rstp, N, regM, supposedVar,
          outOfBound, 0, w);
  delete w;
  if (outOfBound != 0) std::cout << "weight out of bound\n";
  return DataFrame::create(Named("val") = rst, Named("P") = rstp);
}

// first and last are pointint to orderAcc
// return the supposed variance
/*
double divide2sums(std::vector<pointerToDist>::iterator first,
                 std::vector<pointerToDist>::iterator last,
                 std::vector<pointerToDist>::iterator&head,
                 rgMethod&rglinear_add,
                 unsigned&N, unsigned&forDeltaCompare,
                 truncMethod&truncate,
                 double headTrunc, double tailTrunc,

                 std::vector<double>&SD,
                 std::vector<unsigned>&ind,

                 std::vector<std::vector<unsigned> >&gridInd,
                 std::vector<double>&cors,

                 dist&functionValue, unsigned&outOfBound, double&distVariance){

// int first_head=first-head, last_head=last-head;
int indFirst=ind[first-head], indLast=ind[last-head];

if(first+1==last)
{
  dist rst(N);
  formGrid(rst.val.begin(), rst.val.size(),
first->ptr->val.front()+last->ptr->val.front(),
  (first->ptr->val.back()-first->ptr->val.front()+
    (last->ptr->val.back()-last->ptr->val.front()))/(rst.val.size()-1));
  rst.val.back()=first->ptr->val.back()+last->ptr->val.back();

  double cor=0, supposedVar=SD[indLast]*SD[indLast]+SD[indFirst]*SD[indFirst];

  bool pureConv=1;
  if(gridInd.size()!=0)
  {
    pureConv=0;
    for(int i=0,iend=cors.size();i!=iend;++i)
    {
      //if(gridInd[2*i][indFirst]==gridInd[2*i][indLast]&&gridInd[2*i+1][indFirst]==gridInd[2*i+1][indLast])
      if(gridInd[i][indFirst]==gridInd[i][indLast])
      {
        cor=cors[i];
        break;
      }
    }
    supposedVar+=cor*2*SD[indFirst]*SD[indLast];
  }

  distVariance=mixConv(first->ptr->val, first->ptr->P,
SD[indFirst]*SD[indFirst], last->ptr->val, last->ptr->P,
SD[indLast]*SD[indLast], rst.val, rst.P, forDeltaCompare, rglinear_add,
          supposedVar,
          outOfBound, pureConv);
  if(truncate(rst.val, rst.P, headTrunc,
tailTrunc))distVariance=Var(rst.val,rst.P); swapDist(rst,functionValue); return
supposedVar;
}


pointerToDist sum;
sum.range=(first->range+last->range)/2;
std::vector<pointerToDist>::iterator splitPoint=std::lower_bound(first, last+1,
sum);

if(splitPoint->range==sum.range&&splitPoint!=last)
{
  if(splitPoint->range-(splitPoint-1)->range>(splitPoint+1)->range-splitPoint->range)
    ++splitPoint;
}

dist leftBranch, rightBranch, &left=*(first->ptr), &right=*(last->ptr);

double supposedVarLeft=0, supposedVarRight=0, leftDistVar=0, rightDistVar=0;

if(first!=splitPoint-1)
{
  supposedVarLeft=divide2sums(first, splitPoint-1, head,
                         rglinear_add,
                         N, forDeltaCompare,
                         truncate,
                         headTrunc, tailTrunc,
                         SD, ind, gridInd, cors,
                         leftBranch, outOfBound, leftDistVar);
  swapDist(left,leftBranch);
}
else
{
  supposedVarLeft=SD[indFirst]*SD[indFirst];
  leftDistVar=supposedVarLeft;
}

if(splitPoint!=last)
{
  supposedVarRight=divide2sums(splitPoint, last, head,
                          rglinear_add,
                          N, forDeltaCompare,
                          truncate,
                          headTrunc, tailTrunc,
                          SD, ind, gridInd, cors,
                          rightBranch, outOfBound, rightDistVar);
  swapDist(right,rightBranch);
}
else
{
  supposedVarRight=SD[indLast]*SD[indLast];
  rightDistVar=supposedVarRight;
}

double supposedVar=supposedVarLeft+supposedVarRight;

bool pureConv=1;
if(gridInd.size()!=0)
{
  double covar=0;
  pureConv=0;
  for(int i=first-head,iend=splitPoint-head;i!=iend;++i)
  {
    int indI=ind[i];
    for(int j=splitPoint-head,jend=last+1-head;j!=jend;++j)
    {
      int indJ=ind[j], k=0, kend=cors.size();
      for(;k!=kend;++k)
      {
        //if(gridInd[2*k][indI]==gridInd[2*k][indJ]&&gridInd[2*k+1][indI]==gridInd[2*k+1][indJ])
        if(gridInd[k][indI]==gridInd[k][indJ])
        {
          covar+=cors[k]*SD[indI]*SD[indJ];
          break;
        }
      }
      if(k==kend)break;
    }
  }
  supposedVar+=2*covar;
}


dist rst(N);
formGrid(rst.val.begin(), rst.val.size(), left.val.front()+right.val.front(),
        (left.val.back()-left.val.front()+
        (right.val.back()-right.val.front()))/(N-1));
rst.val.back()=left.val.back()+right.val.back();

distVariance=mixConv(left.val, left.P, leftDistVar,
        right.val, right.P, rightDistVar,
        rst.val, rst.P, forDeltaCompare,
        rglinear_add, supposedVar, outOfBound, pureConv);
if(truncate(rst.val, rst.P, headTrunc,
tailTrunc))distVariance=Var(rst.val,rst.P); swapDist(rst,functionValue); return
supposedVar;
}












double divide2sumsWithUnderlyingCor(std::vector<pointerToDist>::iterator first,
                 std::vector<pointerToDist>::iterator last,
                 std::vector<pointerToDist>::iterator&head,
                 rgMethod&rglinear_add,
                 unsigned&N, unsigned&forDeltaCompare,
                 truncMethod&truncate,
                 double headTrunc, double tailTrunc,

                 std::vector<double>&SD,
                 std::vector<unsigned>&ind,

                 std::vector<std::vector<unsigned> >&gridInd,
                 std::vector<double>&cors,

                 double underlyingCor,

                 dist&functionValue, unsigned&outOfBound, double&distVariance){

// int first_head=first-head, last_head=last-head;
int indFirst=ind[first-head], indLast=ind[last-head];

if(first+1==last)
{
  dist rst(N);
  formGrid(rst.val.begin(), rst.val.size(),
first->ptr->val.front()+last->ptr->val.front(),
  (first->ptr->val.back()-first->ptr->val.front()+
    (last->ptr->val.back()-last->ptr->val.front()))/(rst.val.size()-1));
  rst.val.back()=first->ptr->val.back()+last->ptr->val.back();

  double cor=underlyingCor,
    supposedVar=SD[indLast]*SD[indLast]+SD[indFirst]*SD[indFirst];

  // bool pureConv=1;
  // if(gridInd.size()!=0)
  // {
  //   pureConv=0;
    for(int i=0,iend=cors.size();i!=iend;++i)
    {
      //if(gridInd[2*i][indFirst]==gridInd[2*i][indLast]&&gridInd[2*i+1][indFirst]==gridInd[2*i+1][indLast])
      if(gridInd[i][indFirst]==gridInd[i][indLast])
      {
        cor=cors[i];
        break;
      }
    }
    supposedVar+=cor*2*SD[indFirst]*SD[indLast];
  // }

  distVariance=mixConv(first->ptr->val, first->ptr->P,
SD[indFirst]*SD[indFirst], last->ptr->val, last->ptr->P,
SD[indLast]*SD[indLast], rst.val, rst.P, forDeltaCompare, rglinear_add,
          supposedVar,
          outOfBound, 0);
  if(truncate(rst.val, rst.P, headTrunc,
tailTrunc))distVariance=Var(rst.val,rst.P); swapDist(rst,functionValue); return
supposedVar;
}


pointerToDist sum;
sum.range=(first->range+last->range)/2;
std::vector<pointerToDist>::iterator splitPoint=std::lower_bound(first, last+1,
sum);

if(splitPoint->range==sum.range&&splitPoint!=last)
{
  if(splitPoint->range-(splitPoint-1)->range>(splitPoint+1)->range-splitPoint->range)
    ++splitPoint;
}

dist leftBranch, rightBranch, &left=*(first->ptr), &right=*(last->ptr);

double supposedVarLeft=0, supposedVarRight=0, leftDistVar=0, rightDistVar=0;

if(first!=splitPoint-1)
{
  supposedVarLeft=divide2sumsWithUnderlyingCor(first, splitPoint-1, head,
                         rglinear_add,
                         N, forDeltaCompare,
                         truncate,
                         headTrunc, tailTrunc,
                         SD, ind, gridInd, cors, underlyingCor,
                         leftBranch, outOfBound, leftDistVar);
  swapDist(left,leftBranch);
}
else
{
  supposedVarLeft=SD[indFirst]*SD[indFirst];
  leftDistVar=supposedVarLeft;
}

if(splitPoint!=last)
{
  supposedVarRight=divide2sumsWithUnderlyingCor(splitPoint, last, head,
                          rglinear_add,
                          N, forDeltaCompare,
                          truncate,
                          headTrunc, tailTrunc,
                          SD, ind, gridInd, cors, underlyingCor,
                          rightBranch, outOfBound, rightDistVar);
  swapDist(right,rightBranch);
}
else
{
  supposedVarRight=SD[indLast]*SD[indLast];
  rightDistVar=supposedVarRight;
}

double supposedVar=supposedVarLeft+supposedVarRight;

// bool pureConv=1;
// if(gridInd.size()!=0)
// {
  double covar=0,dummyCovar=0;
  // pureConv=0;
  for(int i=first-head,iend=splitPoint-head;i!=iend;++i)
  {
    int indI=ind[i];
    for(int j=splitPoint-head,jend=last+1-head;j!=jend;++j)
    {
      int indJ=ind[j], k=0, kend=cors.size();
      for(;k!=kend;++k)
      {
        //if(gridInd[2*k][indI]==gridInd[2*k][indJ]&&gridInd[2*k+1][indI]==gridInd[2*k+1][indJ])
        if(gridInd[k][indI]==gridInd[k][indJ])
        {
          double tmp=SD[indI]*SD[indJ];
          covar+=cors[k]*tmp;
          dummyCovar+=tmp;
          break;
        }
      }
      if(k==kend)break;
    }
  }

  double sumSDfirstPart=0, sumSDsecondPart=0;
  for(int
i=first-head,iend=splitPoint-head;i!=iend;++i)sumSDfirstPart+=SD[ind[i]];
  for(int
j=splitPoint-head,jend=last+1-head;j!=jend;++j)sumSDsecondPart+=SD[ind[j]];

  supposedVar+=2*(covar+(sumSDfirstPart*sumSDsecondPart-dummyCovar)*underlyingCor);
// }


dist rst(N);
formGrid(rst.val.begin(), rst.val.size(), left.val.front()+right.val.front(),
        (left.val.back()-left.val.front()+
        (right.val.back()-right.val.front()))/(N-1));
rst.val.back()=left.val.back()+right.val.back();

distVariance=mixConv(left.val, left.P, leftDistVar,
        right.val, right.P, rightDistVar,
        rst.val, rst.P, forDeltaCompare,
        rglinear_add, supposedVar, outOfBound, 0);
if(truncate(rst.val, rst.P, headTrunc,
tailTrunc))distVariance=Var(rst.val,rst.P); swapDist(rst,functionValue); return
supposedVar;
}
*/

double divide2sumsNew(std::vector<pointerToDist>::iterator first,
                      std::vector<pointerToDist>::iterator last,
                      std::vector<pointerToDist>::iterator &head,
                      rgMethod &rglinear_add, unsigned &N,
                      unsigned &forDeltaCompare, truncMethod &truncate,
                      double headTrunc, double tailTrunc,

                      std::vector<double> &SD, std::vector<unsigned> &ind,

                      std::vector<std::vector<unsigned> > &gridInd,
                      std::vector<double> &cors,

                      dist &functionValue, unsigned &outOfBound,
                      double &distVariance, wToUseFFT *w) {
  // int first_head=first-head, last_head=last-head;
  int indFirst = ind[first - head], indLast = ind[last - head];

  if (first + 1 == last) {
    dist rst(N);
    formGrid(rst.val.begin(), rst.val.size(),
             first->ptr->val.front() + last->ptr->val.front(),
             (first->ptr->val.back() - first->ptr->val.front() +
              (last->ptr->val.back() - last->ptr->val.front())) /
                 (rst.val.size() - 1));
    rst.val.back() = first->ptr->val.back() + last->ptr->val.back();

    double cor = 0, supposedVar =
                        SD[indLast] * SD[indLast] + SD[indFirst] * SD[indFirst];

    bool pureConv = 1;
    if (gridInd.size() != 0) {
      pureConv = 0;
      for (int i = 0, iend = cors.size(); i != iend; ++i) {
        // if(gridInd[2*i][indFirst]==gridInd[2*i][indLast]&&gridInd[2*i+1][indFirst]==gridInd[2*i+1][indLast])
        if (gridInd[i][indFirst] == gridInd[i][indLast]) {
          cor = cors[i];
          break;
        }
      }
      supposedVar += cor * 2 * SD[indFirst] * SD[indLast];
    }

    distVariance = mixConv(
        first->ptr->val, first->ptr->P, SD[indFirst] * SD[indFirst],
        last->ptr->val, last->ptr->P, SD[indLast] * SD[indLast], rst.val, rst.P,
        forDeltaCompare, rglinear_add, supposedVar, outOfBound, pureConv, w);
    if (truncate(rst.val, rst.P, headTrunc, tailTrunc))
      distVariance = Var(rst.val, rst.P);
    swapDist(rst, functionValue);
    return supposedVar;
  }

  pointerToDist sum;
  sum.range = (first->range + last->range) / 2;
  std::vector<pointerToDist>::iterator splitPoint =
      std::lower_bound(first, last + 1, sum);

  if (splitPoint->range == sum.range && splitPoint != last) {
    if (splitPoint->range - (splitPoint - 1)->range >
        (splitPoint + 1)->range - splitPoint->range)
      ++splitPoint;
  }

  dist leftBranch, rightBranch, &left = *(first->ptr), &right = *(last->ptr);

  double supposedVarLeft = 0, supposedVarRight = 0, leftDistVar = 0,
         rightDistVar = 0;

  if (first != splitPoint - 1) {
    supposedVarLeft =
        divide2sumsNew(first, splitPoint - 1, head, rglinear_add, N,
                       forDeltaCompare, truncate, headTrunc, tailTrunc, SD, ind,
                       gridInd, cors, leftBranch, outOfBound, leftDistVar, w);
    swapDist(left, leftBranch);
  } else {
    supposedVarLeft = SD[indFirst] * SD[indFirst];
    leftDistVar = supposedVarLeft;
  }

  if (splitPoint != last) {
    supposedVarRight =
        divide2sumsNew(splitPoint, last, head, rglinear_add, N, forDeltaCompare,
                       truncate, headTrunc, tailTrunc, SD, ind, gridInd, cors,
                       rightBranch, outOfBound, rightDistVar, w);
    swapDist(right, rightBranch);
  } else {
    supposedVarRight = SD[indLast] * SD[indLast];
    rightDistVar = supposedVarRight;
  }

  double supposedVar = supposedVarLeft + supposedVarRight;

  bool pureConv = 1;
  if (gridInd.size() != 0) {
    double covar = 0;
    pureConv = 0;

    int tmpSize = cors.size();
    double forwardSum[tmpSize], backwardSum[tmpSize], COR[tmpSize];

    COR[tmpSize - 1] = cors[tmpSize - 1];
    for (int i = tmpSize - 2; i >= 0; --i) COR[i] = cors[i] - cors[i + 1];

    for (int i = 0; i != tmpSize; ++i) {
      forwardSum[i] = 0;
      backwardSum[i] = 0;
    }

    int i = splitPoint - head, j = i - 1, k = 0,
        istart = i;  // i is forward, j is backward
    // k is the pointer to the smallest identical grid level
    int iend = last + 1 - head;

    // std::cout<<first-head+1<<" "<<j+1<<" "<<i+1<<" "<<iend<<"----\n";

    while (i != iend) {
      int findk = k;
      for (; findk != tmpSize; ++findk) {
        if (gridInd[findk][ind[i]] == gridInd[findk][ind[j]]) break;
      }
      k = findk;
      if (k == tmpSize) break;
      for (int h = k; h != tmpSize; ++h) forwardSum[h] += SD[ind[i]];
      ++i;
    }
    i = istart;

    int jend = first - head - 1;
    k = 0;
    while (j != jend) {
      int findk = k;
      for (; findk != tmpSize; ++findk) {
        if (gridInd[findk][ind[i]] == gridInd[findk][ind[j]]) break;
      }
      k = findk;
      if (k == tmpSize) break;
      for (int h = k; h != tmpSize; ++h) backwardSum[h] += SD[ind[j]];
      --j;
    }

    for (int i = 0; i != tmpSize; ++i) {
      covar += forwardSum[i] * backwardSum[i] * COR[i];
    }
    // std::cout<<"\n\n";
    supposedVar += 2 * covar;
  }

  dist rst(N);
  formGrid(rst.val.begin(), rst.val.size(),
           left.val.front() + right.val.front(),
           (left.val.back() - left.val.front() +
            (right.val.back() - right.val.front())) /
               (N - 1));
  rst.val.back() = left.val.back() + right.val.back();

  distVariance = mixConv(left.val, left.P, leftDistVar, right.val, right.P,
                         rightDistVar, rst.val, rst.P, forDeltaCompare,
                         rglinear_add, supposedVar, outOfBound, pureConv, w);
  if (truncate(rst.val, rst.P, headTrunc, tailTrunc))
    distVariance = Var(rst.val, rst.P);
  swapDist(rst, functionValue);
  return supposedVar;
}

double closestPairConv(std::vector<pointerToDist> &orderAcc,
                       std::vector<double> &SD, std::vector<unsigned> &ind,
                       std::vector<std::vector<unsigned> > &gridInd,
                       std::vector<double> &cors,
                       // double underlyingCor,
                       char regridMethod, unsigned N, unsigned forDeltaCompare,
                       bool single, double headTrunc,
                       double tailTrunc,  // double cor,
                       dist &functionValue, unsigned &outOfBound,
                       double &mixConvVariance, wToUseFFT *w) {
  rgMethod rglinear_add = nullptr;
  if (regridMethod == 'l')
    rglinear_add = &rglinearAddCpp;
  else
    rglinear_add = &rglr4splitAddCpp;

  truncMethod truncate;
  if (single)
    truncate = &truncateDistSingle;
  else
    truncate = &truncateDist;

  // std::sort(orderAcc.begin(),orderAcc.end());

  for (std::vector<pointerToDist>::iterator i = orderAcc.begin() + 1;
       i != orderAcc.end(); ++i)
    i->range += (i - 1)->range;

  std::vector<pointerToDist>::iterator thehead = orderAcc.begin(),
                                       &head = thehead;

  return divide2sumsNew(orderAcc.begin(), orderAcc.end() - 1, head,
                        rglinear_add, N, forDeltaCompare, truncate, headTrunc,
                        tailTrunc, SD, ind, gridInd, cors, functionValue,
                        outOfBound, mixConvVariance, w);
}

// result is stored in L.front()
// L will be totally polluted
double mixConvOneSample(std::vector<dist> &L, std::vector<unsigned> &ind,
                        std::vector<std::vector<unsigned> > &gridInd,
                        std::vector<double> &cors, char regridMethod,
                        unsigned N, unsigned forDeltaCompare, bool single,
                        double headTrunc, double tailTrunc, dist &functionValue,
                        unsigned &outOfBound, std::vector<double> &SD,
                        wToUseFFT *w) {
  if (L.size() == 1) {
    swapDist(functionValue, L.front());
    return SD[ind.front()] * SD[ind.front()];
  }

  std::vector<pointerToDist> pToL(L.size());
  std::vector<pointerToDist>::iterator pToLi = pToL.begin();
  for (std::vector<dist>::iterator k = L.begin(), kend = L.end(); k != kend;
       ++k, ++pToLi) {
    pToLi->ptr = k;
    pToLi->range = k->val.back() - k->val.front();
  }

  double mixConvVariance;

  return closestPairConv(pToL, SD, ind, gridInd, cors, regridMethod, N,
                         forDeltaCompare, single, headTrunc, tailTrunc,
                         functionValue, outOfBound, mixConvVariance, w);
}

struct parallelInd {
  unsigned start, end;
};

struct experiment : public Worker {
  std::vector<std::vector<unsigned> > &gridInd;
  std::vector<dist> &results;
  std::vector<std::vector<unsigned> > &sampleList;
  std::vector<dist> &fullDistList;
  std::vector<double> &cors;
  std::vector<parallelInd> &theIndex;
  char &regridMethod;
  unsigned &N;
  unsigned &forDeltaCompare;
  bool &single;
  double &headTrunc;
  double &tailTrunc;
  std::vector<unsigned> &outOfBound;
  std::vector<double> &SD, &theoVar;
  wToUseFFT *w;

  // Here is how to initialize class member of references
  experiment(std::vector<std::vector<unsigned> > &gridInd,
             std::vector<dist> &results,
             std::vector<std::vector<unsigned> > &sampleList,
             std::vector<dist> &fullDistList, std::vector<double> &cors,
             std::vector<parallelInd> &theIndex, char &regridMethod,
             unsigned &N, unsigned &forDeltaCompare, bool &single,
             double &headTrunc, double &tailTrunc,
             std::vector<unsigned> &outOfBound, std::vector<double> &SD,
             std::vector<double> &theoVar, wToUseFFT *w)
      : gridInd(gridInd),
        results(results),
        sampleList(sampleList),
        fullDistList(fullDistList),
        cors(cors),
        theIndex(theIndex),
        regridMethod(regridMethod),
        N(N),
        forDeltaCompare(forDeltaCompare),
        single(single),
        headTrunc(headTrunc),
        tailTrunc(tailTrunc),
        outOfBound(outOfBound),
        SD(SD),
        theoVar(theoVar),
        w(w) {}

  void operator()(std::size_t st, std::size_t end) {
    for (std::size_t t = st; t != end; ++t) {
      unsigned i = theIndex[t].start, theEnd = theIndex[t].end;
      for (; i != theEnd; ++i) {
        std::vector<unsigned> &ind =
            sampleList[i];  // ind will be used as rowIndex for gridInd
        std::vector<dist> oneSample(ind.size());

        std::vector<unsigned> indNoZero(ind.size());

        unsigned k = 0;
        for (unsigned j = 0, jend = ind.size(); j != jend; ++j) {
          dist &tmpDist = fullDistList[ind[j] - 1];
          if (tmpDist.val.size() != 1) {
            oneSample[k].val.assign(tmpDist.val.begin(), tmpDist.val.end());
            oneSample[k].P.assign(tmpDist.P.begin(), tmpDist.P.end());
            indNoZero[k] = ind[j] - 1;
            ++k;
          }
        }
        oneSample.resize(k);
        indNoZero.resize(k);
        theoVar[i] =
            mixConvOneSample(oneSample, indNoZero, gridInd, cors, regridMethod,
                             N, forDeltaCompare, single, headTrunc, tailTrunc,
                             results[i], outOfBound[t], SD, w);
      }
    }
  }
};




// after this function, sampleList's elements are sorted
// [[Rcpp::export]]
List mixConvSample(List sampleListR, List fullDistListR, List gridIndex,
                   NumericVector corsR, char regridMethod = 'r',
                   double headTrunc = 1e-100, double tailTrunc = 1e-100,
                   unsigned N = 256, unsigned forDeltaCompare = 256,
                   bool single = 1, int maxCore = 8, bool useFFT = 0) {
  // sampleLis is a list of sorted integer vectors indicating the indexes
  // fullDisLis is a list of data frames, each data frame is a distribution
  // gridInd is a dataframe of integer vectors, each 2 vectors indicate the
  // indexes of cell it should be sorted by the first column to last column the
  // first 2 vectors are for the lowest level cors is the correlation for each
  // level

  std::vector<double> cors(corsR.begin(), corsR.end());

  std::vector<std::vector<unsigned> > gridInd(gridIndex.size() / 2);
  for (int i = 0, iend = gridInd.size(); i != iend; ++i) {
    IntegerVector tmp = gridIndex[2 * i], tmp2 = gridIndex[2 * i + 1];
    gridInd[i].resize(tmp.size());
    unsigned k = 0;
    gridInd[i][0] = k;
    for (unsigned j = 1, jend = tmp.size(); j != jend; ++j) {
      if (tmp[j] != tmp[j - 1] || tmp2[j] != tmp2[j - 1]) ++k;
      gridInd[i][j] = k;
    }
  }

  std::vector<std::vector<unsigned> > sampleList(sampleListR.size());
  for (unsigned i = 0, iend = sampleList.size(); i != iend; ++i) {
    IntegerVector tmp = sampleListR[i];
    sampleList[i].assign(tmp.begin(), tmp.end());
  }

  std::vector<dist> fullDistList(fullDistListR.size());
  for (unsigned i = 0, iend = fullDistList.size(); i != iend; ++i) {
    List tmpDist = fullDistListR[i];
    NumericVector tmpVal = tmpDist[0], tmpP = tmpDist[1];
    if (tmpVal.size() != 2) {
      fullDistList[i].val.assign(tmpVal.begin(), tmpVal.end());
      fullDistList[i].P.assign(tmpP.begin(), tmpP.end());
    } else {
      fullDistList[i].val.resize(3);
      fullDistList[i].val[0] = tmpVal[0];
      fullDistList[i].val[1] = (tmpVal[0] + tmpVal[1]) / 2;
      fullDistList[i].val[2] = tmpVal[1];
      fullDistList[i].P.resize(3);
      fullDistList[i].P[0] = tmpP[0];
      fullDistList[i].P[1] = 0;
      fullDistList[i].P[2] = tmpP[1];
    }
  }

  // std::cout<<"copy done\n";

  // Sleep(30*1000);return List::create();

  std::vector<parallelInd> theIndex(maxCore);
  unsigned interval = sampleList.size() / maxCore;

  for (unsigned i = 0, iend = theIndex.size(); i != iend; ++i) {
    theIndex[i].start = interval * i;
    theIndex[i].end = theIndex[i].start + interval;
  }
  theIndex.back().end = sampleList.size();

  std::vector<unsigned> outOfBound(maxCore, 0);

  std::vector<double> SD(fullDistList.size());

  std::vector<double> theoVar(sampleList.size());

  for (unsigned i = 0, iend = fullDistList.size(); i != iend; ++i) {
    SD[i] = std::sqrt(Var(fullDistList[i].val, fullDistList[i].P));
  }

  // //********

  std::vector<dist> results(sampleList.size());

  wToUseFFT *w = NULL;
  if (useFFT) w = new wToUseFFT(N);

  experiment T(gridInd, results, sampleList, fullDistList, cors, theIndex,
               regridMethod, N, forDeltaCompare, single, headTrunc, tailTrunc,
               outOfBound, SD, theoVar, w);

  parallelFor(0, maxCore, T);
  delete w;

  // std::cout<<"N of out-of-bound weights =
  // "<<std::accumulate(outOfBound.begin(),outOfBound.end(),0)<<"\n";

  List rst(results.size());
  for (int i = 0, iend = results.size(); i != iend; ++i) {
    rst[i] = DataFrame::create(Named("val") = results[i].val,
                               Named("P") = results[i].P);
  }

  return List::create(Named("theoreticalVar") = theoVar,
                      Named("mixConv") = rst);
}


double mixConvOneSampleCovVarSumPerSample(
    std::vector<int> &ind, std::vector<double> &theVar, std::vector<double> &sd,
    std::vector<std::vector<unsigned> > &gridInd, std::vector<double> &cors) {
  double variance = 0;
  for (int i = 0, iend = ind.size(); i != iend; ++i) variance += theVar[ind[i]];

  unsigned corSize = cors.size();
  double cor[corSize];
  cor[corSize - 1] = cors.back();
  for (int i = corSize - 2; i >= 0; --i) cor[i] = cors[i] - cors[i + 1];
  double totalCovar = 0;

  for (int i = 0, iend = gridInd.size(); i != iend; ++i) {
    std::vector<unsigned> &currentInd = gridInd[i];
    int j = 1, jend = ind.size();
    double covarInThisGrid = 0, sumSD = sd[ind[0]], sumVar = theVar[ind[0]];
    unsigned currentIndex, priorIndex = ind[0];
    while (j != jend) {
      currentIndex = ind[j];
      if (currentInd[currentIndex] != currentInd[priorIndex]) {
        covarInThisGrid += (sumSD * sumSD - sumVar);
        sumSD = sd[currentIndex];
        sumVar = theVar[currentIndex];
        priorIndex = currentIndex;
      } else {
        sumSD += sd[currentIndex];
        sumVar += theVar[currentIndex];
      }
      ++j;
    }
    covarInThisGrid += (sumSD * sumSD - sumVar);
    totalCovar += covarInThisGrid * cor[i];
  }
  return variance + totalCovar;
}


struct covarParallel : public Worker {
  std::vector<parallelInd> &theIndex;
  std::vector<std::vector<int> > &sampleList;
  std::vector<double> &theVar;
  std::vector<double> &sd;
  std::vector<std::vector<unsigned> > &gridInd;
  std::vector<double> &cors;
  std::vector<double> &rst;
  // Here is how to initialize class member of references
  covarParallel(std::vector<parallelInd> &theIndex,
                std::vector<std::vector<int> > &sampleList,
                std::vector<double> &theVar, std::vector<double> &sd,
                std::vector<std::vector<unsigned> > &gridInd,
                std::vector<double> &cors, std::vector<double> &rst)
      : theIndex(theIndex),
        sampleList(sampleList),
        theVar(theVar),
        sd(sd),
        gridInd(gridInd),
        cors(cors),
        rst(rst) {}

  void operator()(std::size_t st, std::size_t end) {
    for (std::size_t t = st; t != end; ++t) {
      int i = theIndex[t].start, theEnd = theIndex[t].end;
      for (; i != theEnd; ++i) {
        std::vector<int> &ind = sampleList[i];
        rst[i] =
            mixConvOneSampleCovVarSumPerSample(ind, theVar, sd, gridInd, cors);
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector covarMatSum(List sampleList, NumericVector variances,
                          NumericVector cors, List gridInd, int maxCore = 8) {
  std::vector<double> vars(variances.begin(), variances.end()),
      corrs(cors.begin(), cors.end());
  std::vector<std::vector<int> > SL(sampleList.size());
  for (int i = 0, iend = SL.size(); i != iend; ++i) {
    IntegerVector tmp = sampleList[i];
    IntegerVector tmp2 = tmp - 1;
    SL[i].assign(tmp2.begin(), tmp2.end());
  }

  std::vector<std::vector<unsigned> > GI(gridInd.size() / 2);
  for (int i = 0, iend = GI.size(); i != iend; ++i) {
    IntegerVector tmp = gridInd[2 * i], tmp2 = gridInd[2 * i + 1];
    GI[i].resize(tmp.size());
    unsigned k = 0;
    GI[i][0] = k;
    for (unsigned j = 1, jend = tmp.size(); j != jend; ++j) {
      if (tmp[j] != tmp[j - 1] || tmp2[j] != tmp2[j - 1]) ++k;
      GI[i][j] = k;
    }
  }

  std::vector<parallelInd> theIndex(maxCore);
  int interval = SL.size() / maxCore;
  for (int i = 0, iend = theIndex.size(); i != iend; ++i) {
    theIndex[i].start = interval * i;
    theIndex[i].end = theIndex[i].start + interval;
  }
  theIndex.back().end = sampleList.size();

  std::vector<double> rst(SL.size());
  std::vector<double> sd(vars.size());
  for (unsigned i = 0, iend = sd.size(); i != iend; ++i)
    sd[i] = std::sqrt(vars[i]);
  covarParallel T(theIndex, SL, vars, sd, GI, corrs, rst);

  parallelFor(0, maxCore, T);

  NumericVector result(rst.begin(), rst.end());
  return result;
}

struct levelInfo {
  unsigned nextLevelRowInGridInd;
  double sumSD2, sumSD;
  levelInfo() {
    sumSD2 = 0;
    sumSD = 0;
  }
};

template <typename T>
void sumLevelInfo(T begin, T end, levelInfo &rst) {
  for (; begin != end; ++begin) {
    rst.sumSD += begin->sumSD;
    rst.sumSD2 += (begin->sumSD * begin->sumSD);
  }
}

// [[Rcpp::export]]
NumericVector generateCorMultiplier(NumericVector var, List gridIndR) {
  std::vector<std::vector<unsigned> > gridInd(gridIndR.size() / 2);
  for (int i = 0, iend = gridInd.size(); i != iend; ++i) {
    IntegerVector tmp = gridIndR[2 * i], tmp2 = gridIndR[2 * i + 1];
    gridInd[i].resize(tmp.size());
    unsigned k = 0;
    gridInd[i][0] = k;
    for (unsigned j = 1, jend = tmp.size(); j != jend; ++j) {
      if (tmp[j] != tmp[j - 1] || tmp2[j] != tmp2[j - 1]) ++k;
      gridInd[i][j] = k;
    }
  }

  // find all the 1km grid cell, store sd sum square and sum var
  std::vector<std::vector<levelInfo> > eachLevelInfo(gridInd.size() + 1);
  eachLevelInfo.front().resize(var.size());
  for (unsigned i = 0, iend = var.size(); i != iend; ++i) {
    eachLevelInfo.front()[i].sumSD = std::sqrt(var[i]);
    eachLevelInfo.front()[i].sumSD2 = var[i];
    eachLevelInfo.front()[i].nextLevelRowInGridInd = i;
  }

  for (unsigned i = 1, iend = eachLevelInfo.size(); i != iend; ++i) {
    std::vector<unsigned> &gind = gridInd[i - 1];
    std::vector<levelInfo> &priorLevel = eachLevelInfo[i - 1];

    unsigned siz = 1;
    for (unsigned j = 1, jend = priorLevel.size(); j != jend; ++j) {
      if (gind[priorLevel[j].nextLevelRowInGridInd] !=
          gind[priorLevel[j - 1].nextLevelRowInGridInd])
        ++siz;
    }
    eachLevelInfo[i].resize(siz);
    // gind and eachLevelInfo

    std::vector<levelInfo>::iterator k = eachLevelInfo[i].begin();
    unsigned priorj = 0;
    for (unsigned j = 1, jend = priorLevel.size(); j != jend; ++j) {
      // if(gind[j]!=gind[j-1])
      if (gind[priorLevel[j].nextLevelRowInGridInd] !=
          gind[priorLevel[j - 1].nextLevelRowInGridInd]) {
        sumLevelInfo(priorLevel.begin() + priorj, priorLevel.begin() + j, *k);
        k->nextLevelRowInGridInd = priorLevel[priorj].nextLevelRowInGridInd;
        priorj = j;
        ++k;
      }
    }
    sumLevelInfo(priorLevel.begin() + priorj, priorLevel.end(), *k);
    k->nextLevelRowInGridInd = priorLevel[priorj].nextLevelRowInGridInd;
  }

  NumericVector rst(gridInd.size() + 1);
  rst[0] = std::accumulate(var.begin(), var.end(), 0.0);
  for (unsigned i = 1, iend = eachLevelInfo.size(); i != iend; ++i) {
    std::vector<levelInfo> &currentLevel = eachLevelInfo[i];
    for (unsigned j = 0, jend = currentLevel.size(); j != jend; ++j) {
      rst[i] += (currentLevel[j].sumSD * currentLevel[j].sumSD -
                 currentLevel[j].sumSD2);
    }
  }
  return rst;
}


// {
//
// //--------------------------------------------------------------------------------------------------
// //--------------------------------------------------------------------------------------------------
// //---------------------------mix-convolution, one portfolio. It can be
// applied to multi cor, multi index settings
//
//
// double divide2sumsNewSinglePortfolio(std::vector<pointerToDist>::iterator
// first,
//                  std::vector<pointerToDist>::iterator last,
//                  std::vector<pointerToDist>::iterator&head,
//                  rgMethod&rglinear_add,
//                  unsigned&N, unsigned&forDeltaCompare,
//                  truncMethod&truncate,
//                  double headTrunc, double tailTrunc,
//                  std::vector<double>&SD,
//                  std::vector<std::vector<unsigned> >&gridInd,
//                  std::vector<double>&cors,
//                  dist&functionValue, unsigned&outOfBound,
//                  double&distVariance){
//
//
// int indFirst=first-head, indLast=last-head;
//
// if(first+1==last)
// {
//   dist rst(N);
//   formGrid(rst.val.begin(), rst.val.size(),
//   first->ptr->val.front()+last->ptr->val.front(),
//   (first->ptr->val.back()-first->ptr->val.front()+
//     (last->ptr->val.back()-last->ptr->val.front()))/(rst.val.size()-1));
//   rst.val.back()=first->ptr->val.back()+last->ptr->val.back();
//
//   double cor=0,
//   supposedVar=SD[indLast]*SD[indLast]+SD[indFirst]*SD[indFirst];
//
//   bool pureConv=1;
//   if(gridInd.size()!=0)
//   {
//     pureConv=0;
//     for(int i=0,iend=cors.size();i!=iend;++i)
//     {
//       if(gridInd[i][indFirst]==gridInd[i][indLast])
//       {
//         cor=cors[i];
//         break;
//       }
//     }
//     supposedVar+=cor*2*SD[indFirst]*SD[indLast];
//   }
//
//   distVariance=mixConv(first->ptr->val, first->ptr->P,
//   SD[indFirst]*SD[indFirst],
//           last->ptr->val, last->ptr->P, SD[indLast]*SD[indLast],
//           rst.val, rst.P, forDeltaCompare,
//           rglinear_add,
//           supposedVar,
//           outOfBound, pureConv,w);
//   if(truncate(rst.val, rst.P, headTrunc,
//   tailTrunc))distVariance=Var(rst.val,rst.P); swapDist(rst,functionValue);
//   return supposedVar;
// }
//
//
// pointerToDist sum;
// sum.range=(first->range+last->range)/2;
// std::vector<pointerToDist>::iterator splitPoint=std::lower_bound(first,
// last+1, sum);
//
// if(splitPoint->range==sum.range&&splitPoint!=last)
// {
//   if(splitPoint->range-(splitPoint-1)->range>(splitPoint+1)->range-splitPoint->range)
//     ++splitPoint;
// }
//
// dist leftBranch, rightBranch, &left=*(first->ptr), &right=*(last->ptr);
//
// double supposedVarLeft=0, supposedVarRight=0, leftDistVar=0, rightDistVar=0;
//
// if(first!=splitPoint-1)
// {
//   supposedVarLeft=divide2sumsNewSinglePortfolio(first, splitPoint-1, head,
//                          rglinear_add,
//                          N, forDeltaCompare,
//                          truncate,
//                          headTrunc, tailTrunc,
//                          SD, //ind,
//                          gridInd, cors,
//                          leftBranch, outOfBound, leftDistVar);
//   swapDist(left,leftBranch);
// }
// else
// {
//   supposedVarLeft=SD[indFirst]*SD[indFirst];
//   leftDistVar=supposedVarLeft;
// }
//
// if(splitPoint!=last)
// {
//   supposedVarRight=divide2sumsNewSinglePortfolio(splitPoint, last, head,
//                           rglinear_add,
//                           N, forDeltaCompare,
//                           truncate,
//                           headTrunc, tailTrunc,
//                           SD, //ind,
//                           gridInd, cors,
//                           rightBranch, outOfBound, rightDistVar);
//   swapDist(right,rightBranch);
// }
// else
// {
//   supposedVarRight=SD[indLast]*SD[indLast];
//   rightDistVar=supposedVarRight;
// }
//
// double supposedVar=supposedVarLeft+supposedVarRight;
//
// bool pureConv=1;
// if(gridInd.size()!=0)
// {
//   double covar=0;
//   pureConv=0;
//
//   int tmpSize=cors.size();
//   double forwardSum[tmpSize],backwardSum[tmpSize],COR[tmpSize];
//
//   COR[tmpSize-1]=cors[tmpSize-1];
//   for(int i=tmpSize-2;i>=0;--i)COR[i]=cors[i]-cors[i+1];
//
//   for(int i=0;i!=tmpSize;++i){forwardSum[i]=0;backwardSum[i]=0;}
//
//   int i=splitPoint-head,j=i-1,k=0,istart=i;//i is forward, j is backward
//   // k is the pointer to the smallest identical grid level
//   int iend=last+1-head;
//
//   while(i!=iend)
//   {
//     int findk=k;
//     for(;findk!=tmpSize;++findk)
//     {
//       //if(gridInd[findk][ind[i]]==gridInd[findk][ind[j]])break;
//       if(gridInd[findk][i]==gridInd[findk][j])break;
//     }
//     k=findk;
//     if(k==tmpSize)break;
//     //for(int h=k;h!=tmpSize;++h)forwardSum[h]+=SD[ind[i]];
//     for(int h=k;h!=tmpSize;++h)forwardSum[h]+=SD[i];
//     ++i;
//   }
//   i=istart;
//
//   int jend=first-head-1;
//   k=0;
//   while(j!=jend)
//   {
//     int findk=k;
//     for(;findk!=tmpSize;++findk)
//     {
//       //if(gridInd[findk][ind[i]]==gridInd[findk][ind[j]])break;
//       if(gridInd[findk][i]==gridInd[findk][j])break;
//     }
//     k=findk;
//     if(k==tmpSize)break;
//     //for(int h=k;h!=tmpSize;++h)backwardSum[h]+=SD[ind[j]];
//     for(int h=k;h!=tmpSize;++h)backwardSum[h]+=SD[j];
//     --j;
//   }
//
//   for(int i=0;i!=tmpSize;++i)
//   {
//     covar+=forwardSum[i]*backwardSum[i]*COR[i];
//   }
//   supposedVar+=2*covar;
// }
//
//
// dist rst(N);
// formGrid(rst.val.begin(), rst.val.size(), left.val.front()+right.val.front(),
//         (left.val.back()-left.val.front()+
//         (right.val.back()-right.val.front()))/(N-1));
// rst.val.back()=left.val.back()+right.val.back();
//
// distVariance=mixConv(left.val, left.P, leftDistVar,
//         right.val, right.P, rightDistVar,
//         rst.val, rst.P, forDeltaCompare,
//         rglinear_add, supposedVar, outOfBound, pureConv);
// if(truncate(rst.val, rst.P, headTrunc,
// tailTrunc))distVariance=Var(rst.val,rst.P); swapDist(rst,functionValue);
// return supposedVar;
// }
//
//
//
//
//
//
//
//
// double closestPairConvSinglePortfolio(std::vector<pointerToDist>&orderAcc,
//                     std::vector<double>&SD,
//                     std::vector<std::vector<unsigned> >&gridInd,
//                     std::vector<double>&cors,
//                     char regridMethod,
//                     unsigned N, unsigned forDeltaCompare, bool single,
//                     double headTrunc, double tailTrunc, //double cor,
//                     dist&functionValue, unsigned&outOfBound,
//                     double&mixConvVariance){
//
// rgMethod rglinear_add = nullptr;
// if(regridMethod=='l')rglinear_add=&rglinearAddCpp;
// else rglinear_add=&rglr4splitAddCpp;
//
// truncMethod truncate;
// if(single)truncate=&truncateDistSingle;
// else truncate=&truncateDist;
//
// //std::sort(orderAcc.begin(),orderAcc.end());
//
// for(std::vector<pointerToDist>::iterator
// i=orderAcc.begin()+1;i!=orderAcc.end();++i)
//   i->range+=(i-1)->range;
//
// std::vector<pointerToDist>::iterator thehead=orderAcc.begin(), &head=thehead;
//
// /*
// if(underlyingCor<1e-10)return
// divide2sums(orderAcc.begin(),orderAcc.end()-1,head,
//                    rglinear_add,N,forDeltaCompare,
//                    truncate, headTrunc, tailTrunc,
//                    SD, ind, gridInd, cors,
//                    functionValue,outOfBound,mixConvVariance);
//
// return divide2sumsWithUnderlyingCor(orderAcc.begin(),orderAcc.end()-1,head,
//                    rglinear_add,N,forDeltaCompare,
//                    truncate, headTrunc, tailTrunc,
//                    SD, ind, gridInd, cors, underlyingCor,
//                    functionValue,outOfBound,mixConvVariance);
// */
//
//
// return divide2sumsNewSinglePortfolio(orderAcc.begin(),orderAcc.end()-1,head,
//                    rglinear_add,N,forDeltaCompare,
//                    truncate, headTrunc, tailTrunc,
//                    SD, //ind,
//                    gridInd, cors,
//                    functionValue,outOfBound,mixConvVariance);
// }
//
//
//
//
//
//
// double mixConvSinglePortfolio(std::vector<dist>&L,
// //std::vector<unsigned>&ind,
//                       std::vector<std::vector<unsigned> >&gridInd,
//                       std::vector<double>&cors,
//                       char regridMethod, unsigned N, unsigned
//                       forDeltaCompare, bool single, double headTrunc, double
//                       tailTrunc, dist&functionValue, unsigned&outOfBound,
//                       std::vector<double>&SD){
//
//
// //std::vector<dist>L=Li;
//
// if(L.size()==1)
// {
//   swapDist(functionValue,L.front());
//   //return SD[ind.front()]*SD[ind.front()];
//   return SD[0]*SD[0];
// }
//
//
// std::vector<pointerToDist>pToL(L.size());
// std::vector<pointerToDist>::iterator pToLi=pToL.begin();
// for(std::vector<dist>::iterator k=L.begin(),kend=L.end();k!=kend;++k,++pToLi)
// {
//   pToLi->ptr=k;
//   pToLi->range=k->val.back()-k->val.front();
// }
//
// double mixConvVariance;
// return closestPairConvSinglePortfolio(pToL, SD, //ind,
//                                       gridInd, cors,
//                        regridMethod, N, forDeltaCompare,
//                                    single, headTrunc, tailTrunc,
//                                    functionValue, outOfBound,
//                                    mixConvVariance);
//
// }
//
//
//
//
//
//
// struct singlePortfolioMulCor: public Worker{
// std::vector<std::vector<unsigned> >&gridInd;
// std::vector<dist>&results;
// std::vector<dist>&fullDistList;
// std::vector<std::vector<double> >&cors;
// std::vector<parallelInd>&theIndex;
// char&regridMethod;
// unsigned&N;
// unsigned&forDeltaCompare;
// bool&single;
// double&headTrunc;
// double&tailTrunc;
// std::vector<unsigned>&outOfBound;
// std::vector<double>&SD, &theoVar;
//
// //Here is how to initialize class member of references
// singlePortfolioMulCor(std::vector<std::vector<unsigned> >&gridInd,
//             std::vector<dist>&results,
//             std::vector<dist>&fullDistList,
//             std::vector<std::vector<double> >&cors,
//             std::vector<parallelInd>&theIndex,
//             char&regridMethod,
//             unsigned&N,
//             unsigned&forDeltaCompare,
//             bool&single,
//             double&headTrunc,
//             double&tailTrunc,
//             std::vector<unsigned>&outOfBound,
//             std::vector<double>&SD,
//             std::vector<double>&theoVar):
//               gridInd(gridInd), results(results),
//               fullDistList(fullDistList),
//               cors(cors),theIndex(theIndex),
//               regridMethod(regridMethod),N(N),forDeltaCompare(forDeltaCompare),
//               single(single),headTrunc(headTrunc),tailTrunc(tailTrunc),outOfBound(outOfBound),
//               SD(SD), theoVar(theoVar){}
//
// void operator()(std::size_t st, std::size_t end){
// for(std::size_t t=st;t!=end;++t)
// {
//   unsigned i=theIndex[t].start, theEnd=theIndex[t].end;
//   for(;i!=theEnd;++i)
//   {
//       std::vector<dist>copyfullDistList=fullDistList;
//       theoVar[i]=mixConvSinglePortfolio(fullDistList, gridInd, cors[i],
//       regridMethod, N, forDeltaCompare, single,
//                       headTrunc, tailTrunc, results[i], outOfBound[t], SD);
//   }
// }
// }
// };
//
//
//
//
//
//
// // [[Rcpp::export]]
// List mixConvSinglePortfolioMultiCors(List fullDistListR, List gridIndex, List
// corsR,
//                    char regridMethod='r',
//                    double headTrunc=1e-100, double tailTrunc=1e-100,
//                    unsigned N=256, unsigned forDeltaCompare=256, bool
//                    single=1, int maxCore=8){
// // sampleLis is a list of sorted integer vectors indicating the indexes
// // fullDisLis is a list of data frames, each data frame is a distribution
// // gridInd is a dataframe of integer vectors, each 2 vectors indicate the
// indexes of cell
// // it should be sorted by the first column to last column
// // the first 2 vectors are for the lowest level
// // cors is the correlation for each level
//
// std::vector<std::vector<double> >cors(corsR.size());
// for(int i=0,iend=corsR.size();i!=iend;++i)
// {
//   NumericVector tmp=corsR[i];
//   cors[i].assign(tmp.begin(),tmp.end());
// }
//
// //std::vector<double>cors(corsR.begin(),corsR.end());
//
//
// std::vector<std::vector<unsigned> >gridInd(gridIndex.size()/2);
// for(int i=0,iend=gridInd.size();i!=iend;++i)
// {
//   IntegerVector tmp=gridIndex[2*i], tmp2=gridIndex[2*i+1];
//   gridInd[i].resize(tmp.size());
//   unsigned k=0;
//   gridInd[i][0]=k;
//   for(unsigned j=1,jend=tmp.size();j!=jend;++j)
//   {
//     if(tmp[j]!=tmp[j-1]||tmp2[j]!=tmp2[j-1])++k;
//     gridInd[i][j]=k;
//   }
// }
//
//
// std::vector<dist>fullDistList(fullDistListR.size());
// unsigned tmpSize=fullDistList.size();
// unsigned ind[tmpSize];
// unsigned k=0;
// for(unsigned i=0,iend=fullDistList.size();i!=iend;++i)
// {
//   List tmpDist=fullDistListR[i];
//   NumericVector tmpVal=tmpDist[0], tmpP=tmpDist[1];
//   if(tmpVal.size()==1)continue;
//   ind[k]=i;
//   fullDistList[k].val.assign(tmpVal.begin(),tmpVal.end());
//   fullDistList[k].P.assign(tmpP.begin(),tmpP.end());
//   ++k;
// }
// fullDistList.resize(k);
//
//
// for(unsigned i=0,iend=gridInd.size();i!=iend;++i)
// {
//   std::vector<unsigned>&currentGridInd=gridInd[i];
//   for(unsigned j=0;j!=k;++j)currentGridInd[j]=currentGridInd[ind[j]];
//   currentGridInd.resize(k);
// }
//
//
// std::cout<<"copy done\n";
//
// std::vector<parallelInd>theIndex(maxCore);
// unsigned interval=cors.size()/maxCore;
//
// for(unsigned i=0,iend=theIndex.size();i!=iend;++i)
// {
//   theIndex[i].start=interval*i;
//   theIndex[i].end=theIndex[i].start+interval;
// }
// theIndex.back().end=cors.size();
//
// std::vector<unsigned>outOfBound(maxCore,0);
//
// std::vector<double>SD(fullDistList.size());
//
// std::vector<double>theoVar(cors.size());
//
// for(unsigned i=0,iend=fullDistList.size();i!=iend;++i)
// {
//   SD[i]=std::sqrt(Var(fullDistList[i].val,fullDistList[i].P));
// }
//
// // //********
// std::vector<dist>results(cors.size());
//
// // mixConvSinglePortfolio(fullDistList, gridInd, cors[0], regridMethod, N,
// forDeltaCompare, single,
// //                  headTrunc, tailTrunc, results[0], outOfBound[0], SD);
// //
// // return List::create();
//
// singlePortfolioMulCor T(gridInd, results, fullDistList, cors, theIndex,
//             regridMethod, N, forDeltaCompare, single, headTrunc, tailTrunc,
//             outOfBound, SD, theoVar);
//
// parallelFor(0, maxCore, T);
//
// std::cout<<"N of out-of-bound weights =
// "<<std::accumulate(outOfBound.begin(),outOfBound.end(),0)<<"\n";
//
// List rst(results.size());
// for(int i=0,iend=results.size();i!=iend;++i)
// {
//   rst[i]=DataFrame::create(Named("val")=results[i].val,
//   Named("P")=results[i].P);
// }
//
// return List::create(Named("theoreticalVar")=theoVar,Named("mixConv")=rst);
// }
//
//
//
//
// }

// struct dist{
// std::vector<losstype>val;
// std::vector<probtype>P;
// dist(unsigned N)
// {
//   val.resize(N);
//   P.resize(N,0);
// }
// dist(){};
// void getGrid(losstype start, losstype end)
// {
//   val.front()=start;
//   losstype d=(end-start)/(val.size()-1);
//   for(std::vector<losstype>::iterator i=val.begin()+1;i!=val.end();++i)
//     *i=*(i-1)+d;
//   val.back()=end;
//   std::fill(P.begin(), P.end(), 0);
// }
// };


// rgMethod="l" or "r"
void comoDistListcpp(dist &result, vec<dist> &distList, char rgMethod,
                     unsigned maxPoint) {
  dist rstPrior(maxPoint * 2);
  dist rst(maxPoint * 2);

  comoVec(distList.front().val, distList.front().P, distList[1].val, distList[1].P,
       rstPrior.val, rstPrior.P);

  bool isRegriddedInTheEnd = 0;

  for (unsigned i = 2, iend = distList.size(); i < iend; ++i) {
    // losstype theMax=rstPrior.val.back()+distList[i].val.back();

    if (rstPrior.val.size() + distList[i].val.size() > maxPoint) {
      dist tmpRst(rstPrior.val.size() + distList[i].val.size());
      comoVec(distList[i].val, distList[i].P, rstPrior.val, rstPrior.P, tmpRst.val,
           tmpRst.P);
      rstPrior.val.resize(maxPoint);
      rstPrior.P.resize(maxPoint);
      rstPrior.getGrid(tmpRst.val.front(), tmpRst.val.back());
      std::fill(rstPrior.P.begin(), rstPrior.P.end(), 0);
      if (rgMethod == 'l') {
        rglinearAddCpp(tmpRst.val.begin(), tmpRst.P.begin(), tmpRst.val.size(),
                       rstPrior.val.begin(), rstPrior.P.begin(),
                       rstPrior.val.size());
      } else {
        rglr4splitAddCpp(tmpRst.val.begin(), tmpRst.P.begin(),
                         tmpRst.val.size(), rstPrior.val.begin(),
                         rstPrior.P.begin(), rstPrior.val.size());
      }
      isRegriddedInTheEnd = 1;
    } else {
      rst.P.resize(maxPoint);
      rst.val.resize(maxPoint);
      comoVec(distList[i].val, distList[i].P, rstPrior.val, rstPrior.P, rst.val,
           rst.P);
      swapDist(rstPrior, rst);
      isRegriddedInTheEnd = 0;
    }

    // Rcout<<theMax<<" "<<rstPrior.val.back()<<"\n";
  }

  if (!isRegriddedInTheEnd) {
    unsigned tmpint = rstPrior.val.size();
    rst.val.resize(std::min(tmpint, maxPoint));
    rst.P.resize(rst.val.size());
    std::fill(rst.P.begin(), rst.P.end(), 0);
    rst.getGrid(rstPrior.val.front(), rstPrior.val.back());

    if (rgMethod == 'l') {
      rglinearAddCpp(rstPrior.val.begin(), rstPrior.P.begin(),
                     rstPrior.val.size(), rst.val.begin(), rst.P.begin(),
                     rst.val.size());
    } else {
      rglr4splitAddCpp(rstPrior.val.begin(), rstPrior.P.begin(),
                       rstPrior.val.size(), rst.val.begin(), rst.P.begin(),
                       rst.val.size());
    }
    swapDist(rst, rstPrior);
  }

  swapDist(result, rstPrior);
}







//'
//' Comonotonize a list of distributions
//' 
//' Comonotonize a list of PMFs sequentially.
//' 
//' @inheritParams sequentialConvAll
//' 
//' @param regridMethod 
//' Regriding method:
//' 
//'   \code{"lr"}: linear regrid. Default.
//' 
//'   \code{"r4"}: four-point regrid.
//' 
//' @inheritParams convSplitAtom4productsIrregular
//' 
//' @details  See \code{\link{como}}.
//' 
//' @inherit rglr return
//' 
//' @example inst/examples/comoDistList.R
//'
// [[Rcpp::export]]
DataFrame comoDistList(
    List lisOfDists, 
    String regridMethod = "lr",
    int N = 256) 
{
  
  
  checkPMFlistIntegrity(lisOfDists);
  
  
  vec<dist> distList(lisOfDists.size());
  for (unsigned i = 0, iend = distList.size(); i != iend; ++i) {
    List X = lisOfDists[i];
    NumericVector val = X[0], P = X[1];
    distList[i].val.assign(val.begin(), val.end());
    distList[i].P.assign(P.begin(), P.end());
  }

  dist result;
  char rgd = regridMethod == "lr" ? 'l' : 'r';
  comoDistListcpp(result, distList, rgd, N);

  return DataFrame::create(Named("val") = result.val, Named("P") = result.P);
}




// // [[Rcpp::export]]
// DataFrame mixConvSim(List distList, List gridInd, NumericVector cors,
//                      NumericVector distVarsR=-1,
//                      std::string convMethod="nearestNeighbor",
//                      double headTruncInConv=1e-100, double
//                      TailTruncInConv=1e-100, double headTruncInEnd=1e-10,
//                      double tailTruncInEnd=1e-10, char regridMethod='l', int
//                      N=256, int forDeltaCompare=256, bool single=1, int
//                      maxCore=7){
//
//
// vec<dist>DL(distList.size());
// {
//   for(unsigned i=0,iend=DL.size();i!=iend;++i)
//   {
//     List tmp=distList[i];
//     NumericVector val=tmp[0], P=tmp[1];
//     DL[i].val.assign(val.begin(),val.end());
//     DL[i].P.assign(P.begin(),P.end());
//   }
// }
//
//
// vec<losstype>distVars;
// losstype*vars=&*distVarsR.begin();
// if(distVarsR[0]==-1)
// {
//   distVars.resize(DL.size());
//   for(unsigned i=0,iend=distVars.size();i!=iend;++i)
//   {
//     distVars[i]=Var(DL[i].val, DL[i].P);
//   }
//   vars=&distVars.front();
// }
//
//
//
//
//
//
//
//
//
// }
//
//
//
//
//

// spit out recursive nearest neighbor convolutional order
losstype rnnConvOrder(losstype *accWidthStart, indtype currentStart,
                      indtype currentEnd, indtype *&pairIndex,
                      indtype *layerInd, indtype layerN,
                      indtype **spatialCellBounds, losstype *layerCor,
                      losstype *accS, losstype *sd, losstype *&theoVar,
                      unsigned char &hierachy, unsigned char *hier)
// The number of exposures is 'N'.

// 'layerN' is the number of hierarchies of spatial correlation model.

// 'layerInd' is an integer vector of size 'N'*'layerN'. Each of 'layerN'
// elements are spatial cell indexes for an exposure.

// 'spatialCellBounds' is a size-'layerN' vector of pointers.
// 'spatialCellBounds[i]' represents an integer vector.
// 'spatialCellBounds[i][k]' is the index of the 1st exposure located in the
// k_th spatial cell of the i_th hierachy. 'spatialCellBounds[i]' is of size k+1

// 'spatialCellBoundsSizes' represents a size-'layerN' integer vector. The i_th
// element is the size of vector spatialCellBounds[i]

{
  if (currentEnd - currentStart == 1) {
    hierachy = -1 + 1;
    return sd[currentStart] * sd[currentStart];
    // return the variance
  }

  // losstype*rightBegin=NULL;
  indtype rightBegin = 0;
  {
    // losstype mid=(accWidth[siz-1]+*(accWidth-1))/2;
    losstype mid =
        (accWidthStart[currentEnd - 1] + accWidthStart[(int)currentStart - 1]) /
        2;

    rightBegin = std::lower_bound(accWidthStart + currentStart,
                                  accWidthStart + currentEnd, mid) -
                 accWidthStart;
    if (rightBegin != currentEnd - 1 &&
        accWidthStart[rightBegin] + accWidthStart[rightBegin - 1] < 2 * mid) {
      ++rightBegin;
    }

    *pairIndex = currentEnd - 1;
    --pairIndex;
    *pairIndex = rightBegin - 1;
    --pairIndex;
  }

  // for the right part
  unsigned char rightHierachy;
  losstype rightVar = rnnConvOrder(
      accWidthStart, rightBegin, currentEnd, pairIndex, layerInd, layerN,
      spatialCellBounds, layerCor, accS, sd, theoVar, rightHierachy, hier - 1);

  // for the left part
  unsigned char leftHierachy;
  losstype leftVar = rnnConvOrder(
      accWidthStart, currentStart, rightBegin, pairIndex, layerInd, layerN,
      spatialCellBounds, layerCor, accS, sd, theoVar, leftHierachy, hier - 1);

  *hier = std::max(leftHierachy, rightHierachy);

  // for variance aggregation
  indtype *layerIndex = layerInd + rightBegin * layerN;
  losstype covar = 0;
  for (indtype h = 0; h != layerN; ++h) {
    indtype cellIndex = layerIndex[h];
    int lb = spatialCellBounds[h][cellIndex];
    int lbEnd = spatialCellBounds[h][cellIndex + 1];
    int mid = rightBegin;
    covar += (accS[lbEnd - 1] - accS[mid - 1]) *
             (accS[mid - 1] - accS[lb - 1]) * layerCor[h];
  }

  losstype theoreticalVar = leftVar + 2 * covar + rightVar;
  *theoVar = theoreticalVar;
  --theoVar;
  ++hierachy;
  return theoreticalVar;
}

// // theoVar is a vector of size N-1. theoVar points to the last element
// void extractOrder(losstype*accSupportWidth, indtype N, losstype*sd,
// indtype*layerInd, losstype*theoVar,
//                   indtype layerN, indtype**spatialCellBounds,
//                   losstype*layerCor, vec<vec<indtype> >&hierOrder)
// {
//   vec<losstype>accS(N+1,0);
//   std::partial_sum(sd, sd+N, &accS[1]);
//   vec<unsigned char>nodeHier(N-1);
//   unsigned char*hierback=&nodeHier.back();
//   // hier is going to be fulfilled
//   vec<indtype>pairIndex((N-1)*2);
//   indtype*pairIndexBack=&pairIndex.back();
//   unsigned char hierachy;
//
//
//   rnnConvOrder(accSupportWidth, 0, N, pairIndexBack, layerInd, layerN,
//   spatialCellBounds,
//                layerCor, &accS[0], sd, theoVar, hierachy, hierback);
//
//
//   hierOrder.resize(hierachy);
//   {
//     vec<indtype>eachHierSize(hierachy,0);
//     {
//       for(indtype i=0,iend=nodeHier.size();i!=iend;++i)
//       {
//         ++eachHierSize[nodeHier[i]];
//       }
//     }
//     for(indtype i=0,iend=hierOrder.size();i!=iend;++i)
//     {
//       hierOrder[i].reserve(eachHierSize[i]*2);
//     }
//     for(indtype i=0,iend=nodeHier.size();i!=iend;++i)
//     {
//       hierOrder[nodeHier[i]].push_back(pairIndex[i*2]);
//       hierOrder[nodeHier[i]].push_back(pairIndex[i*2+1]);
//     }
//   }
// }
//
//
//
//
// // [[Rcpp::export]]
// IntegerVector giveRnnConvOrder(NumericVector supportWidths)
// {
//   IntegerVector rst((supportWidths.size()-1)*2);
//   std::vector<losstype>sw(supportWidths.size()+1,0);
//   std::partial_sum(supportWidths.begin(),supportWidths.end(),sw.begin()+1);
//   losstype*swBegin=&sw[0]+1;
//   unsigned*rstLast=(unsigned*)&*rst.end()-1;
//   rnnConvOrder(swBegin, sw.size()-1, rstLast, swBegin);
//   return rst+1;
// }






























//'
//' Distribution mean
//' 
//' Compute the mean of a PMF.
//' 
//' @inheritParams rglr
//' 
//' @return Mean of the PMF
//'
// [[Rcpp::export]]
double Mean(List X) 
{
  NumericVector val = X[0], P = X[1];
  return std::inner_product(val.begin(), val.end(), P.begin(), 0.0);
}


//'
//' Distribution variance
//' 
//' Compute the variance of a PMF.
//' 
//' @inheritParams rglr
//' 
//' @return Variance of the PMF
//'
// [[Rcpp::export]]
double Var(List X) 
{
  NumericVector val = X[0], P = X[1];
  double m = 0, m2 = 0;
  for(int i = 0, iend = val.size(); i < iend; ++i)
  {
    double tmp = val[i] * P[i];
    m += tmp;
    m2 += tmp * val[i];
  }
  return m2 - m * m;
}












































