#include <Rcpp.h>
using namespace Rcpp;
#include "../h/stratifiedSample.hpp"
#include "../h/miniPCG.hpp"


// [[Rcpp::export]]
IntegerVector testStratifiedSampling(int popuSize, int sampleSize, 
                                     bool replace = false, int seed = 42) 
{
  IntegerVector rst(sampleSize);
  MiniPcg32 rng(seed);
  StratifiedSampling()(popuSize, sampleSize, &rst[0], replace, rng);
  return rst;
}

