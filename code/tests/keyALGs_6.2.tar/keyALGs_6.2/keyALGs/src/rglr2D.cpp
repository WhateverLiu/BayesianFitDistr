// #pragma once


#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
DataFrame rglr2Dcpp(DataFrame dist, NumericVector ngx, NumericVector ngy){
// dist has been sorted by the grid cell index by ngx by ngy
NumericVector x=dist[0], y=dist[1], p=dist[2], ngp(ngx.size()*ngy.size(),0.0);

int i=1,j=1,iend=ngx.size(),jend=ngy.size();
for(int k=0,kend=x.size();k!=kend;++k)
{
  if(x[k]>ngx[i])j=1;
  while(i!=iend&&x[k]>ngx[i])++i;
  if(y[k]<ngy[j-1])j=1;
  while(j!=jend&&y[k]>ngy[j])++j;

  double left=x[k]-ngx[i-1], right=ngx[i]-x[k],
         down=y[k]-ngy[j-1], up=ngy[j]-y[k], tmp=(left+right)*(down+up)/p[k];

  int downLeftInd=(i-1)*jend+(j-1), upLeftInd=downLeftInd+1,
    downRightInd=downLeftInd+jend, upRightInd=downRightInd+1;
  ngp[downLeftInd]+=right*up/tmp;
  ngp[upLeftInd]+=right*down/tmp;
  ngp[downRightInd]+=left*up/tmp;
  ngp[upRightInd]+=left*down/tmp;
}

NumericVector gridX(iend*jend), gridY(iend*jend);
for(int k=0,kend=ngx.size();k!=kend;++k)
  std::fill(gridX.begin()+k*jend, gridX.begin()+k*jend+jend, ngx[k]);

for(int k=0,kend=ngx.size();k!=kend;++k)
  std::copy(ngy.begin(),ngy.end(),gridY.begin()+k*jend);

return DataFrame::create(Named("x")=gridX,Named("y")=gridY,Named("P")=ngp);
}









