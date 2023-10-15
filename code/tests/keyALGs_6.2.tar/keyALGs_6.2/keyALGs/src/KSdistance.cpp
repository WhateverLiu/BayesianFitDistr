// #pragma once

#include <Rcpp.h>
using namespace Rcpp;


double distance(std::vector<double>&x, std::vector<double>&y, String ms)
{
std::vector<double> z(x.size());
unsigned i=0, end=z.size();
double tmp, tmp1;
if(ms=="Cheb")
{
  for(;i!=end;++i)z[i]=std::abs(x[i]-y[i]);
  return *std::max_element(z.begin(), z.end());
}
if(ms=="Euc")
{
  for(;i!=end;++i)z[i]=pow(x[i]-y[i],2);
  return sqrt(std::accumulate(z.begin(), z.end(), 0.0)/end);
}
if(ms=="CB")
{
  for(;i!=end;++i)z[i]=std::abs(x[i]-y[i]);
  return std::accumulate(z.begin(), z.end(), 0.0)/end;
}
if(ms=="Can")
{
  for(;i!=end;++i)
  {
    tmp=x[i]+y[i];
    if(tmp==0.0)z[i]=0;
    else z[i]=std::abs(x[i]-y[i])/tmp;
  }
  return std::accumulate(z.begin(), z.end(), 0.0);
}
if(ms=="B")
{
  for(;i!=end;++i)z[i]=sqrt(x[i]*y[i]);
  tmp=std::accumulate(z.begin(), z.end(), 0.0);
  if(tmp<=0)return 0.0;
  return -log(tmp);
}
if(ms=="H")
{
  for(;i!=end;++i)z[i]=sqrt(x[i]*y[i]);
  return 2*(1-sqrt((std::accumulate(z.begin(), z.end(), 0.0))));
}
if(ms=="J")
{
  for(;i!=end;++i)
  {
    if(x[i]<=0.0||y[i]<=0.0)z[i]=0;
    else z[i]=(x[i]-y[i])*log(x[i]/y[i]);
  }
  return std::accumulate(z.begin(), z.end(), 0.0);
}
if(ms=="Top")
{
  for(;i!=end;++i)
  {
    if(x[i]<=0.0)tmp=0; else tmp=x[i]*log(2*x[i]/(x[i]+y[i]));
    if(y[i]<=0.0)tmp1=0; else tmp1=y[i]*log(2*y[i]/(x[i]+y[i]));
    z[i]=tmp1+tmp;
  }
  return std::accumulate(z.begin(), z.end(), 0.0);
}
if(ms=="Tani")
{
  std::vector<double>z1(z.size());
  for(;i!=end;++i)
  {
    z[i]=std::max(x[i], y[i]);
    z1[i]=std::min(x[i], y[i]);
  }
  return 1-std::accumulate(z1.begin(), z1.end(), 0.0)/std::accumulate(z.begin(), z.end(), 0.0);
}
std::cout<<"Hey, you forgot to specify the metrics.\n";
return 0;
}


// [[Rcpp::export]]
double vecD(NumericVector x, NumericVector y, String ms="Cheb")
{
  if(x.size()!=y.size()){std::cout<<"length not the same!!\n";return 0;}
  std::vector<double>xx=as<std::vector<double> >(x),yy=as<std::vector<double> >(y);
  return distance(xx,yy,ms);
}






