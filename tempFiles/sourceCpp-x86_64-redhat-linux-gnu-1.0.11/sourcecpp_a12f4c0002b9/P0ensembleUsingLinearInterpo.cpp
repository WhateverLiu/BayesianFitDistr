// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "h/charlieThreadPool2.hpp"
#include "h/monoLinearInterpo.hpp"
#include "h/miniPCG.hpp"
// #include "h/stratifiedSample.hpp"
#define vec std::vector


template<typename ing, typename num>
struct SampleDataMakeP0
{
  vec<ing> ind;
  vec<num> mdr, dr;
  void sampleData(num *fullMDR, num *fullCDR, ing Nfull,
                  MiniPcg32 &rng, ing sampleSize)
  {
    ind.resize(Nfull);
    std::iota(ind.begin(), ind.end(), 0);
    std::shuffle(ind.begin(), ind.end(), rng);
    std::sort(ind.begin(), ind.begin() + sampleSize);
    mdr.resize(sampleSize);
    dr.resize(sampleSize);
    for (ing i = 0, iend = mdr.size(); i < iend; ++i)
    {
      mdr[i] = fullMDR[ind[i]];
      dr[i] = fullCDR[ind[i]];
    }
  }
  
  
  // Run after sampleData()
  vec<ing> zeroCumsum;
  vec<num> mdrCumsum, Mmdr, p0;
  void makeP0(num windowSizeRatio, num slidingSpeedRatio, num eps) 
  {
    // vec<int> zeroCumsum(dr.size() + 1);
    int windowSize = std::round(windowSizeRatio * mdr.size());
    int slidingSpeed = std::round(slidingSpeedRatio * mdr.size());
    
    
    zeroCumsum.resize(dr.size() + 1);
    zeroCumsum[0] = 0;
    int *cs = zeroCumsum.data() + 1;
    for (int i = 0, iend = dr.size(); i < iend; ++i)
      cs[i] = cs[i - 1] + (dr[i] < eps);
    
    
    // vec<double> mdrCumsum(mdr.size() + 1);
    mdrCumsum.resize(mdr.size() + 1);
    mdrCumsum[0] = 0;
    std::partial_sum(mdr.begin(), mdr.end(), mdrCumsum.begin() + 1);
    double *mdrCumsumPtr = mdrCumsum.data() + 1;
    
    
    int Npoint = mdr.size() / slidingSpeed + 3;
    // vec<double> Mmdr(Npoint), p0(Npoint);
    Mmdr.resize(Npoint);
    p0.resize(Npoint);
    Mmdr[0] = 0;
    p0[0] = 1;
    Mmdr.resize(1);
    p0.resize(1);
    
    
    int mdrSize = mdr.size();
    int begin = 0;
    while (true)
    {
      int end = begin + windowSize;
      if (end >= mdrSize)
      {
        int tmp = std::max(-1, mdrSize - 1 - windowSize); // Position before window.
        int tmp2 = mdrSize - 1 - tmp; // Actual window size.
        double mmdr = (mdrCumsumPtr[mdrSize - 1] - mdrCumsumPtr[tmp]) / tmp2;
        double tmpP0 = (cs[mdrSize - 1] - cs[tmp]) / (tmp2 + 0.0);
        Mmdr.emplace_back(mmdr);
        p0.emplace_back(tmpP0);
        break;
      }
      double mmdr = (mdrCumsumPtr[end - 1] - mdrCumsumPtr[begin - 1]) / windowSize;
      double tmpP0 = (cs[end - 1] - cs[begin - 1]) / (windowSize + 0.0);
      Mmdr.emplace_back(mmdr);
      p0.emplace_back(tmpP0);
      begin += slidingSpeed;
    }
    
    
    Mmdr.emplace_back(1);
    p0.emplace_back(0);
  }
  
  
  // Build the function.
  MonoLinearInterpo<ing, num, num> mlt;
  // By now, the only containers that matter are Mmdr, p0.
  // available vec<num> are: mdrCumsum, mdr, dr.
  // available vec<ing> are: ind, zeroCumsum
  void makeFunAndPredict(num *targetMDRs, num *targetMDRsEnd, num *rst)
  {
    mlt.xsubseq.swap(mdr);
    mlt.ysubseq.swap(dr);
    mlt.indx.swap(ind);
    mlt.lms.C.swap(zeroCumsum);
    mlt.lms.subseq.swap(mdrCumsum);
    
    
    mlt.template reset<3> (
        &*Mmdr.begin(), &*Mmdr.end(), &*p0.begin());
    mlt(targetMDRs, targetMDRsEnd, rst);
    
    
    mlt.xsubseq.swap(mdr);
    mlt.ysubseq.swap(dr);
    mlt.indx.swap(ind);
    mlt.lms.C.swap(zeroCumsum);
    mlt.lms.subseq.swap(mdrCumsum);
  }
    
  
};


// [[Rcpp::export]]
List p0ensembleLDSlinear(NumericVector rawMdr, 
                         NumericVector dr, // NumericVector P0,
                         int windowSize,
                         int slidingSpeed,
                         double eps,
                         NumericVector targetMDRs,
                         int Nsample,
                         double sampleSize,
                         int seed,
                         int maxCore) 
{
  CharlieThreadPool cp(std::move(maxCore));
  vec<SampleDataMakeP0<int, double> > smpV(maxCore);
  vec<MiniPcg32> rngV(maxCore);
  for (int i = 0; i < maxCore; ++i) rngV[i].seed(seed + i);
  double windowSizeRatio = double(windowSize) / dr.size();
  double slidingSpeedRatio = double(slidingSpeed) / dr.size();
  vec<double*> rst(Nsample);
  List result(Nsample);
  for (int i = 0; i < Nsample; ++i)
  {
    NumericVector v(targetMDRs.size());
    result[i] = v;
    rst[i] = &v[0];
  }
  auto f = [&](std::size_t i, std::size_t t)->bool
  {
    smpV[t].sampleData(&rawMdr[0], &dr[0], dr.size(), rngV[t], sampleSize);
    smpV[t].makeP0(windowSizeRatio, slidingSpeedRatio, eps);
    smpV[t].makeFunAndPredict(&targetMDRs[0], &*targetMDRs.end(), rst[i]);
    return false;
  };
  cp.parFor(0, Nsample, f, 1);
  return result;
}  












// mdr should have been shorted. mdr.size() == dr.size()
vec<vec<double> > makeP0(
    NumericVector mdr, NumericVector dr, 
    int windowSize, int slidingSpeed, double eps)
{
  vec<int> zeroCumsum(dr.size() + 1);
  zeroCumsum[0] = 0;
  int *cs = zeroCumsum.data() + 1;
  for (int i = 0, iend = dr.size(); i < iend; ++i)
    cs[i] = cs[i - 1] + (dr[i] < eps);

  
  vec<double> mdrCumsum(mdr.size() + 1);
  mdrCumsum[0] = 0;
  std::partial_sum(mdr.begin(), mdr.end(), mdrCumsum.begin() + 1);
  double *mdrCumsumPtr = mdrCumsum.data() + 1;
  
  
  int Npoint = mdr.size() / slidingSpeed + 1;
  vec<double> Mmdr(Npoint), p0(Npoint);
  Mmdr.resize(0);
  p0.resize(0);
  
  
  int mdrSize = mdr.size();
  int begin = 0;
  while (true)
  {
    int end = begin + windowSize;
    if (end >= mdrSize)
    {
      int tmp = std::max(-1, mdrSize - 1 - windowSize); // Position before window.
      int tmp2 = mdrSize - 1 - tmp; // Actual window size.
      double mmdr = (mdrCumsumPtr[mdrSize - 1] - mdrCumsumPtr[tmp]) / tmp2;
      double tmpP0 = (cs[mdrSize - 1] - cs[tmp]) / (tmp2 + 0.0);
      Mmdr.emplace_back(mmdr);
      p0.emplace_back(tmpP0);
      break;
    }
    double mmdr = (mdrCumsumPtr[end - 1] - mdrCumsumPtr[begin - 1]) / windowSize;
    double tmpP0 = (cs[end - 1] - cs[begin - 1]) / (windowSize + 0.0);
    Mmdr.emplace_back(mmdr);
    p0.emplace_back(tmpP0);
    begin += slidingSpeed;
  }
  
  
  vec<vec<double> > rst(2);
  rst[0].swap(Mmdr);
  rst[1].swap(p0);
  return rst;
}


// [[Rcpp::export]]
List makeP0test(
  NumericVector mdr, NumericVector dr, 
  int windowSize, int slidingSpeed, double eps)
{
  auto rst = makeP0(
    mdr, dr, 
    windowSize, slidingSpeed, eps);
  return List::create(Named("MDR") = rst[0],
                      Named("P0") = rst[1]);
}


// mdr has been sorted.
// [[Rcpp::export]]
List p0ensembleLDSlinearOld(NumericVector rawMdr, 
                         NumericVector dr, // NumericVector P0,
                         int windowSize,
                         int slidingSpeed,
                         double eps,
                         NumericVector targetMDRs,
                         int Nsample,
                         double sampleSizeRatio,
                         int seed,
                         int maxCore) 
{
  
  
  List p0models(Nsample);
  vec<double*> samples(Nsample); // Pointers to the result.
  for (int i =0; i < Nsample; ++i)
  {
    NumericVector v(targetMDRs.size());
    samples[i] = &v[0];
    p0models[i] = v;
  }
  
  
  auto mdrAndP0 = makeP0(
    rawMdr, dr, windowSize, slidingSpeed, eps);
  vec<double> &mdr = mdrAndP0[0], &P0 = mdrAndP0[1];
  int popuSize = mdr.size();
  int sampleSize = std::max<int> (1, std::round(sampleSizeRatio * popuSize));
  
  
  CharlieThreadPool cp(std::move(maxCore));
  vec<MiniPcg32> rngs(maxCore);
  vec<vec<int> > inds(maxCore, vec<int>(0));
  for (int i = 0; i < maxCore; ++i) rngs[i].seed(seed + i);
  vec<MonoLinearInterpo<int, double, double> > monoIntpV(maxCore);
  vec<vec<double> > sampledMDRs(maxCore);
  vec<vec<double> > sampledP0s(maxCore);
  
  
  auto f = [&](std::size_t i, std::size_t t)->bool
  {
    auto &rng = rngs[t];
    auto &ind = inds[t];
    auto &sampledMDR = sampledMDRs[t];
    auto &sampledP0 = sampledP0s[t];
    auto rst = samples[i];
    
    
    if (ind.size() == 0) // Store the sampled indices.
    {
      ind.resize(popuSize);
      std::iota(ind.begin(), ind.end(), 0);
    }
    std::shuffle(ind.begin(), ind.end(), rng);
    std::sort(ind.begin(), ind.begin() + sampleSize);
    
    
    sampledMDR.resize(sampleSize + 2);
    sampledP0.resize(sampleSize + 2);
    sampledMDR[0] = 0;
    sampledP0[0] = 1;
    
    
    // Rcout << "ind.size() = " << ind.size() << "\n";
    // Rcout << "sampleSize = " << sampleSize << "\n";
    
    
    for (int i = 0; i < sampleSize; ++i)
    {
      auto k = ind[i];
      sampledMDR[i + 1] = mdr[k];
      sampledP0[i + 1] = P0[k];
    }
    sampledMDR.back() = 1;
    sampledP0.back() = 0;
    
    
    auto &monoIntp = monoIntpV[t];
    monoIntp.template reset<3> (
        &*sampledMDR.begin(), &*sampledMDR.end(), &*sampledP0.begin());
    monoIntp(&*targetMDRs.begin(), &*targetMDRs.end(), rst);
    
    
    return false;
  };
  
  
  cp.parFor(0, Nsample, f, 1);
  NumericVector mn(targetMDRs.size(), 0.0);
  for (int i = 0; i < Nsample; ++i)
  {
    auto p = samples[i];
    for (int j = 0, jend = targetMDRs.size(); j < jend; ++j) mn[j] += p[j];
  }
  for (int j = 0, jend = targetMDRs.size(); j < jend; ++j) mn[j] /= Nsample;
  
  
  // Rcout << mn[6] << "\n";
  return List::create(Named("p0models") = p0models,
                      Named("P0") = mn);
  // return List::create(Named("p0models") = mn);
}




#undef vec

#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// p0ensembleLDSlinear
List p0ensembleLDSlinear(NumericVector rawMdr, NumericVector dr, int windowSize, int slidingSpeed, double eps, NumericVector targetMDRs, int Nsample, double sampleSize, int seed, int maxCore);
RcppExport SEXP sourceCpp_32_p0ensembleLDSlinear(SEXP rawMdrSEXP, SEXP drSEXP, SEXP windowSizeSEXP, SEXP slidingSpeedSEXP, SEXP epsSEXP, SEXP targetMDRsSEXP, SEXP NsampleSEXP, SEXP sampleSizeSEXP, SEXP seedSEXP, SEXP maxCoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rawMdr(rawMdrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dr(drSEXP);
    Rcpp::traits::input_parameter< int >::type windowSize(windowSizeSEXP);
    Rcpp::traits::input_parameter< int >::type slidingSpeed(slidingSpeedSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type targetMDRs(targetMDRsSEXP);
    Rcpp::traits::input_parameter< int >::type Nsample(NsampleSEXP);
    Rcpp::traits::input_parameter< double >::type sampleSize(sampleSizeSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type maxCore(maxCoreSEXP);
    rcpp_result_gen = Rcpp::wrap(p0ensembleLDSlinear(rawMdr, dr, windowSize, slidingSpeed, eps, targetMDRs, Nsample, sampleSize, seed, maxCore));
    return rcpp_result_gen;
END_RCPP
}
// makeP0test
List makeP0test(NumericVector mdr, NumericVector dr, int windowSize, int slidingSpeed, double eps);
RcppExport SEXP sourceCpp_32_makeP0test(SEXP mdrSEXP, SEXP drSEXP, SEXP windowSizeSEXP, SEXP slidingSpeedSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mdr(mdrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dr(drSEXP);
    Rcpp::traits::input_parameter< int >::type windowSize(windowSizeSEXP);
    Rcpp::traits::input_parameter< int >::type slidingSpeed(slidingSpeedSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(makeP0test(mdr, dr, windowSize, slidingSpeed, eps));
    return rcpp_result_gen;
END_RCPP
}
// p0ensembleLDSlinearOld
List p0ensembleLDSlinearOld(NumericVector rawMdr, NumericVector dr, int windowSize, int slidingSpeed, double eps, NumericVector targetMDRs, int Nsample, double sampleSizeRatio, int seed, int maxCore);
RcppExport SEXP sourceCpp_32_p0ensembleLDSlinearOld(SEXP rawMdrSEXP, SEXP drSEXP, SEXP windowSizeSEXP, SEXP slidingSpeedSEXP, SEXP epsSEXP, SEXP targetMDRsSEXP, SEXP NsampleSEXP, SEXP sampleSizeRatioSEXP, SEXP seedSEXP, SEXP maxCoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rawMdr(rawMdrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dr(drSEXP);
    Rcpp::traits::input_parameter< int >::type windowSize(windowSizeSEXP);
    Rcpp::traits::input_parameter< int >::type slidingSpeed(slidingSpeedSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type targetMDRs(targetMDRsSEXP);
    Rcpp::traits::input_parameter< int >::type Nsample(NsampleSEXP);
    Rcpp::traits::input_parameter< double >::type sampleSizeRatio(sampleSizeRatioSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type maxCore(maxCoreSEXP);
    rcpp_result_gen = Rcpp::wrap(p0ensembleLDSlinearOld(rawMdr, dr, windowSize, slidingSpeed, eps, targetMDRs, Nsample, sampleSizeRatio, seed, maxCore));
    return rcpp_result_gen;
END_RCPP
}
