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



#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// testStratifiedSampling
IntegerVector testStratifiedSampling(int popuSize, int sampleSize, bool replace, int seed);
RcppExport SEXP sourceCpp_28_testStratifiedSampling(SEXP popuSizeSEXP, SEXP sampleSizeSEXP, SEXP replaceSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type popuSize(popuSizeSEXP);
    Rcpp::traits::input_parameter< int >::type sampleSize(sampleSizeSEXP);
    Rcpp::traits::input_parameter< bool >::type replace(replaceSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(testStratifiedSampling(popuSize, sampleSize, replace, seed));
    return rcpp_result_gen;
END_RCPP
}
