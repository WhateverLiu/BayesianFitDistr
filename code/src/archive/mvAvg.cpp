// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector mvAvg(NumericVector x, int start = 1, 
                    int windowSize = 1, int speed = 1) 
{
  start = std::min<int> (std::max<int> (0, start - 1), x.size() - 1);
  double *v = &x[0] + start;
  int size = &*x.end() - v;
  if (windowSize >= size)
    return std::accumulate(v, v + size, 0.0) / size;
  int rstSize = size / speed, r = size % speed;
  if (r == 0) rstSize -= 1;
  rstSize -= r == 0;
  
  
  
}


