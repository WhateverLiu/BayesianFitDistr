// #pragma once


// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;


// #include "h/checkPMFinput.hpp"
extern void checkPMFintegrity(Rcpp::List X, bool checkReverseOrder = false);
extern void checkPMFlistIntegrity(Rcpp::List distlist, bool checkReverseOrder = false);



void rglinearPtr(double *x_st, double *xp_st, unsigned sizex,
                 double *ng_st, double *ngp_st, unsigned ngsize)
{
  unsigned j = 0;
  *ngp_st = 0.0;
  if ( *(ng_st + ngsize-1) < *(x_st + sizex - 1) )
    *(ng_st + ngsize - 1) = *(x_st + sizex - 1);
  double tmp, interval;
  for (unsigned i = 1; i < ngsize; ++i)
  {
    *(ngp_st + i) = 0.0;
    interval = *(ng_st + i) - *(ng_st + i - 1);
    while ( j < sizex and *(x_st + j) <= *(ng_st + i))
    {
      tmp = *(xp_st + j) / interval * ( *(ng_st + i) - *(x_st + j) );
      *(ngp_st + i - 1) = *(ngp_st + i - 1) + tmp;
      *(ngp_st + i) = *(ngp_st + i) + *(xp_st + j) - tmp;
      ++j;
    }
  }
}



void rglinear(std::vector<double>::iterator x_st, std::vector<double>::iterator xp_st,
              unsigned sizex, std::vector<double>::iterator ng_st,
              std::vector<double>::iterator ngp_st, unsigned ngsize = 64)
  // warning!: this function will check ng_st.back() and *(x_st+sizex-1) immediately
  // if ng_st.back() is less, it will set it to *(x_st+sizex-1)!
{
  rglinearPtr(&*x_st,  &*xp_st, sizex, &*ng_st, &*ngp_st, ngsize);
}



void rglr4splitAddCppInThisFile(std::vector<double>::iterator x_st,
                                std::vector<double>::iterator xp_st,
                                unsigned sizex, std::vector<double>::iterator ng_st,
                                std::vector<double>::iterator ngp_st_final,
                                unsigned ngsize)
{
  // warning!: this function will check ng_st.back() and *(x_st+sizex-1) immediately
  // if ng_st.back() is less, it will set it to *(x_st+sizex-1)!

  if(*(ng_st+ngsize-1)<*(x_st+sizex-1))*(ng_st+ngsize-1)=*(x_st+sizex-1);
  std::vector<double>ngpTmp(ngsize,0);
  std::vector<double>::iterator ngp_st=ngpTmp.begin();

  unsigned j=0;

  std::vector<double>::iterator endOfNgp=ngp_st+ngsize-1, endOfNg=ng_st+ngsize-1;
  std::vector<double>::iterator ng_stReserve=ng_st;
  double x0, x1, x2, x3, x, p, p0, p1, p2, p3;
  bool indicator=0;
  // if the distribution being regridded has size < 5 or the resulting distribution
  if(ngsize<5||sizex<5)indicator=1;

  if(!indicator)
  {
    for(unsigned i=1;i!=ngsize;++i)
    {
      x1=*(ng_st+i-1);
      x2=*(ng_st+i);
      while(j < sizex and *(x_st+j)<=*(ng_st+i))
      {
        x3=*endOfNg;
        x0=*ng_st;
        x=*(x_st+j);
        p=*(xp_st+j);
        if(std::abs(x2-x)<1e-20)
        {
          *(ngp_st+i)+=p;
        }
        else if(std::abs(x1-x)<1e-20)
        {
          *(ngp_st+i-1)+=p;
        }
        else
        {
          p1=p*(x3-x)*(x0-x)*(x2-x)/((x3-x1)*(x0-x1)*(x2-x)+(x3-x2)*(x0-x2)*(x-x1));
          p2=(x-x1)*p1/(x2-x);
          p0=(p*(x3-x)-(x3-x1)*p1-(x3-x2)*p2)/(x3-x0);
          p3=p-p0-p1-p2;
          *(ngp_st+i-1)+=p1;
          *(ngp_st+i)+=p2;
          *endOfNgp+=p3;
          *ngp_st+=p0;
        }
        ++j;
      }
    }
    while(1)
    {
      if(*ngp_st<0)
      {
        x=*ng_st;
        x1=*(ng_st+1);
        x2=*(ng_st+2);
        x3=*endOfNg;
        p=*ngp_st;
        p1=(x-x3)*(x-x2)/((x1-x3)*(x1-x2))*p;
        p2=(x-x1)*(x-x3)/((x2-x3)*(x2-x1))*p;
        p3=p-p1-p2;
        *(ngp_st+1)+=p1;
        *(ngp_st+2)+=p2;
        *endOfNgp+=p3;
        ++ng_st;
        *ngp_st=0;
        ++ngp_st;
      }
      if(*ngp_st>=0&&*endOfNgp>=0)break;
      if(endOfNg-ng_st<3)
      {
        indicator=1;
        break;
      }
      if(*endOfNgp<0)
      {
        x=*endOfNg;
        x1=*(endOfNg-1);
        x2=*(endOfNg-2);
        x3=*ng_st;
        p=*endOfNgp;
        p1=(x-x3)*(x-x2)/((x1-x3)*(x1-x2))*p;
        p2=(x-x1)*(x-x3)/((x2-x3)*(x2-x1))*p;
        p3=p-p1-p2;
        *(endOfNgp-1)+=p1;
        *(endOfNgp-2)+=p2;
        *ngp_st+=p3;
        --endOfNg;
        *endOfNgp=0;
        --endOfNgp;
      }
      if(*ngp_st>=0&&*endOfNgp>=0)break;
      if(endOfNg-ng_st<3)
      {
        indicator=1;
        break;
      }
    }
  }

  if(indicator)
  {
    unsigned j=0;
    ng_st=ng_stReserve;
    ngp_st=ngp_st_final;
    double tmp, interval;
    for(unsigned i=1;i!=ngsize;++i)
    {
      interval=*(ng_st+i)-*(ng_st+i-1);
      while(j < sizex and *(x_st+j)<=*(ng_st+i))
      {
        tmp=*(xp_st+j)/interval*(*(ng_st+i)-*(x_st+j));
        *(ngp_st+i-1)+=tmp;
        *(ngp_st+i)+=*(xp_st+j)-tmp;
        ++j;
      }
    }
  }
  else
  {
    ++endOfNgp;
    for(ngp_st_final+=ngp_st-ngpTmp.begin();ngp_st!=endOfNgp;++ngp_st,++ngp_st_final)
      *ngp_st_final+=*ngp_st;
  }
}




//'
//' Linear regriding
//'
//' Linear regrid a PMF
//'
//' @param X A 2-column data frame \code{(support, probability)} as the PMF.
//' Elements in \code{support} must be sorted and have no duplicates.
//' \code{probability} should be nonnegative.
//'
//' @param ngd A sorted numeric vector as the new grid: \code{ngd[1] <= support[1]}
//' AND \code{ngd[length(ngd)] >= support[length(support)]}.
//'
//' @details See \href{https://www.mdpi.com/2227-9091/7/2/54/htm}{Direct and
//' Hierarchical Models for Aggregating Spatially Dependent Catastrophe Risks}.
//'
//' @return A 2-column data frame as the regrided PMF. Column names are
//' \code{"val"} and \code{"P"}.
//'
//' @example inst/examples/rglr.R
//'
// [[Rcpp::export]]
DataFrame rglr(DataFrame X, NumericVector ngd)
{
  checkPMFintegrity(X);
  NumericVector x = X[0], p = X[1], newgridP(ngd.size());
  rglinearPtr(&*x.begin(), &*p.begin(), x.size(), &*ngd.begin(),
              &*newgridP.begin(), ngd.size());
  return DataFrame::create(Named("val") = ngd, Named("P") = newgridP);
}




unsigned findGrid2(double Min, double Max, std::vector<double>::iterator rst_st,
                   unsigned rst_size)
{
  if(Max-Min<1e-12)return 0;
  double interval=(Max-Min)/(rst_size-1);
  std::vector<double>::iterator i=rst_st, end=i+rst_size-1;
  *i=Min;
  ++i;
  for(;i!=end;++i)*i=*(i-1)+interval;
  *i=Max;
  return 1;
}




//'
//' Brute-force convolution
//'
//' Brute-force convolve two PMFs where each PMF's support can be fully irregular.
//' 
//' @inheritParams rglr
//' 
//' @param Y  Same as \code{X}.
//' 
//' @param eps  Numeric threshold to define equality between double-precision 
//' floats: if \code{|a - b| < eps}, then \code{a} and \code{b} are treated as 
//' equal. Default \code{1e-10}. 
//' 
//' @inherit rglr return
//' 
//' @example inst/examples/convBt.R
//' 
// [[Rcpp::export]]
DataFrame convBt(DataFrame X, DataFrame Y, double eps = 1e-10)
{
  checkPMFintegrity(X);
  checkPMFintegrity(Y);
  
  
  NumericVector x = X[0], xp = X[1], y = Y[0], yp = Y[1];
  typedef std::pair<double, double> vp;
  std::vector<vp> distr(x.size() * y.size());
  int k = 0;
  for(int i = 0, iend = x.size(); i != iend; ++i)
  {
    for(int j = 0, jend = y.size(); j != jend; ++j)
    {
      distr[k].first = x[i] + y[j];
      distr[k].second = xp[i] * yp[j];
      ++k;
    }
  }


  std::sort(distr.begin(), distr.end(), [](const vp &x, const vp &y)->bool
  {
    return x.first < y.first;
  });


  // Aggregate probabilities associated with equal value points.
  if (true)
  {
    int i = 0, j = 1, jend = distr.size();
    while (j < jend)
    {
      if (distr[j].first - distr[i].first < eps) // Support values are equal.
      {
        distr[i].second += distr[j].second;
        ++j;
      }
      else
      {
        ++i;
        distr[i] = distr[j];
        ++j;
      }
    }
    distr.resize(i + 1);
  }


  NumericVector rst(distr.size()), rstp(distr.size());
  for (int i = 0, iend = rst.size(); i < iend; ++i)
  {
    rst [i] = distr[i].first;
    rstp[i] = distr[i].second;
  }
  return DataFrame::create(Named("val") = rst, Named("P") = rstp);
}




// the inverse of linear regridding
// [[Rcpp::export]]
DataFrame inverseRglr(DataFrame X, NumericVector ngd, bool normalize = true)
{
  NumericVector val = X[0], P = X[1];
  int distSize = val.size();
  if (val[0] < ngd[0] or val[distSize - 1] > ngd[ngd.size() - 1])
  {
    Rcpp::Rcout << "New grid is narrower than the original. Wrong!\n";
    return DataFrame::create();
  }


  NumericVector ngdP(ngd.size());


  int t = 0, ngdSize = ngd.size(), i = 1;
  while (val[i] < ngd[0]) ++i;
  while (ngd[t] < val[0]) ++t;
  for (int iend = val.size(); i < iend; ++i)
  {
    double interval = val[i] - val[i - 1];
    while (t < ngdSize and ngd[t] <= val[i])
    {
      double w = (ngd[t] - val[i - 1]) / interval;
      ngdP[t] = (1 - w) * P[i - 1] + w * P[i];
      ++t;
    }
  }
  if(normalize) return DataFrame::create(Named("val") = ngd,
     Named("P") = ngdP / std::accumulate(ngdP.begin(), ngdP.end(), 0.0));


  return DataFrame::create(Named("val") = ngd, Named("P") = ngdP);
}





































