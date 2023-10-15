// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "h/trbCharlie.hpp"
#include "h/charlieThreadPool2.hpp"
#include "h/regrid.h"
#include "h/upscaleBoundedMean.h"


#define vec std::vector


// Fine-tuned discretization for fulfilling QA requirements. ftdfqa
struct Ftdfqa
{
  CharlieThreadPool cp;
  std::vector<TrB> trbvec;
  std::vector<Regrid> regridervec;
  
  
  void reset(int &&maxCore) 
  {
    cp.reset(std::move(maxCore)); 
    trbvec.resize(cp.maxCore); 
    regridervec.resize(cp.maxCore);
    maxIterMaxPointCntr.resize(cp.maxCore * 2, 0);
  }
  
  
  Ftdfqa(){ reset(10000); }
  Ftdfqa(int maxCore) {  reset(std::move(maxCore)); }
  Ftdfqa(int &&maxCore) {  reset(std::move(maxCore)); }
  
  
  template <int RIBlib>
  void fineDiscretizeOneTrB(
      double mx, double p0, double a, double b, double c, double d, 
      double *P, int size, TrB &trb)
  {
    trb.reset(a, b, c, d);
    double delta = mx / (size - 1);
    P[0] = p0;
    double halfdelta = delta / 2;
    double priorCdf = trb.cdf<RIBlib> (halfdelta);
    for (int i = 1, iend = size - 1; i < iend; ++i)
    {
      double currentCdf = trb.cdf<RIBlib> (i * delta + halfdelta);
      P[i] = currentCdf - priorCdf;
      priorCdf = currentCdf;
    } 
    P[size - 1] = trb.cdf<RIBlib> ((size - 1) * delta - halfdelta, 0); // tail prob.
    double nr = (1.0 - p0) / std::accumulate(P + 1, P + size, 0.0);
    for (int i = 1; i < size; ++i) P[i] *= nr;
  } 
  
  
  // Impose MDR by linearly scaling the maxes.
  // Return the new size.
  // int maxIterationsForImposingMean, maxPointsRemoved;
  std::vector<int> maxIterMaxPointCntr;
  template <bool downScaleSupport>
  double imposeMDR(
      double *val, double *P, int &size, double MDR, double xmax,
      std::size_t t // which thread this function is being called in.
    )
  {
    int iter = 0;
    double r = upscaleBoundedMean(val, P, size, MDR, xmax, 1e-12, iter, 100);
    int oldsize = size;
    
    
    // maxIterMaxPointCntr.resize(cp.maxCore, 0);
    int *maxIterationsForImposingMean = &maxIterMaxPointCntr[0];
    int *maxPointsRemoved = &maxIterMaxPointCntr[cp.maxCore];
    
    
    maxIterationsForImposingMean[t] = std::max(
      maxIterationsForImposingMean[t], iter);
    
    
    // Notice that, here the final regriding to the 42/64-point will re-impose
    //   the max. But that method could make the second last probability
    //   unreasonably high.
    if (r <= 1)
    {
      if (downScaleSupport)
      {
        for (int i = 0; i < size; ++i) val[i] *= r;  
      }
      else // Adjust P0 to fulfill the requirement.
      {
        P[0] = 1.0 - (1.0 - P[0]) * r;
        for (int i = 1; i < size; ++i) P[i] *= r;
      }
    }
    else
    {
      double sump = 0;
      for (int i = 0; i < size; ++i)
      {
        sump += P[i];
        val[i] = std::min(val[i] * r, xmax);
        if (val[i] >= xmax)
        {
          size = i + 1;
          P[i] += 1.0 - sump;
          break;
        }
      }
    }
    
    
    // if (cp.maxCore <= 1)
    // {
    //   maxPointsRemoved = std::max(maxPointsRemoved, size - oldsize);
    // }
    maxPointsRemoved[t] = std::max(maxPointsRemoved[t], size - oldsize);
    
    
    if constexpr (false)
    {
      double mean = std::inner_product(val, val + size, P, 0.0);
      Rcout << mean << ", " << MDR << ", " << r << "\n";
    }
    
    
    return r;
  }
  
  
  std::vector<double> cntr;
  // TrBtable is a matrix that has 7 rows: MDR, max, p0, a, b, c, d
  template <int regridMethod, int RIBlib, bool downScaleSupport>
  void discretize(double *TrBtable, int ncol, int fineDiscretizeSize, 
                  // double *distTable, 
                  double **distTables, // distTables[j] will be changed!
                  // int supportSize
                  int *supportSizes,
                  int NsupportSizes // Number of different support sizes.
                    )
    // Here, ncol or Ndist is 18999, i.e. the # of nondegenerate PMFs, for example,
    // distTable[] is a space of (supportSize + 2) * 19001
    // The first and last distributions are enclosures of degenerate PMFs.
    // We first fine-discretize TrBs using, e.g., 2000 points. The result
    // is then regrided onto the final support.
  {
    int &Ndist = ncol;
    int maxSuptSize = *std::max_element(supportSizes, supportSizes + NsupportSizes);
    
    
    // Start the table with a degenerate distribution.
    for (int j =0; j < NsupportSizes; ++j)
    {
      std::fill(distTables[j], distTables[j] + supportSizes[j] + 2, 0);
      distTables[j][2] = 1;
      distTables[j] += supportSizes[j] + 2;
    }
    
    
    cntr.resize(cp.maxCore * (fineDiscretizeSize * 2 + maxSuptSize) );
    auto f = [&](std::size_t i, std::size_t t)->bool
    {
      // double *val = cntr.data() + t * (fineDiscretizeSize * 2 + supportSize);
      double *val = cntr.data() + t * (fineDiscretizeSize * 2 + maxSuptSize);
      double *P = val + fineDiscretizeSize;
      double *support = P + fineDiscretizeSize;
      double *param = TrBtable + i * 7;
      makeRegularGrid(0.0, param[1], val, fineDiscretizeSize);
      fineDiscretizeOneTrB<RIBlib> (
          param[1], param[2], param[3], param[4], param[5], param[6], // val, 
          P, fineDiscretizeSize, trbvec[t]);
      
      
      if constexpr (false)
      {
        double mean = std::inner_product(val, val + fineDiscretizeSize, P, 0.0);
        Rcout << mean << ", " << param[0] << "\n";
      }
      
      
      int sizeAfterMDRimposed = fineDiscretizeSize; // sizeAfterMDRimposed is referenced in imposeMDR.
      imposeMDR<downScaleSupport> (val, P, sizeAfterMDRimposed, param[0], param[1], t);
      
      
      double valmax = param[1];
      for (int j = 0; j < NsupportSizes; ++j)
      {
        makeRegularGrid(0.0, valmax, support, supportSizes[j]);  
        double *entry = distTables[j] + i * (supportSizes[j] + 2);
        regridervec[t].operator()<
          regridMethod, int, double, double, double, double, false, true>(
              val, P, sizeAfterMDRimposed, 
              support, entry + 2, supportSizes[j]);
        entry[0] = param[0];
        entry[1] = valmax;
      }
      
      
      return false;
    };
    auto &m = cp.maxCore;
    int grainSize = Ndist / (m*m*m) + 1;
    
    
    std::fill(maxIterMaxPointCntr.begin(), maxIterMaxPointCntr.end(), 0);
    cp.parFor(0, Ndist, f, grainSize);
    
    
    // Check if the maxes are increasing. If not, regrid the current PMF onto
    //   the support of the previous PMF.
    for (int j = 0; j < NsupportSizes; ++j)
    {
      auto distTable = distTables[j];
      auto supportSize = supportSizes[j];
      for (int64_t i = 1, iend = Ndist - 1; i < iend; ++i)
      {
        double &oldmax = distTable[(i - 1) * (supportSize + 2) + 1];
        double &nowmax = distTable[i * (supportSize + 2) + 1];
        if (oldmax <= nowmax) continue;
        double *P = &nowmax + 1;
        double *oldval = cntr.data();
        double *nowval = oldval + supportSize;
        double *oldP = nowval + supportSize;
        makeRegularGrid(0, nowmax, oldval, supportSize);
        makeRegularGrid(0, oldmax, nowval, supportSize);
        std::copy(P, P + supportSize, oldP);
        regridervec[0].operator()<
          regridMethod, int, double, double, double, double, false, true>(
              oldval, oldP, supportSize, nowval, P, supportSize);
        nowmax = oldmax;  
      }  
    }
    
    
    // End the table with degenerate distribution. 
    for (int j = 0; j < NsupportSizes; ++j)
    {
      auto distTable = distTables[j];
      auto supportSize = supportSizes[j];
      distTable += int64_t(supportSize + 2) * Ndist;
      std::fill(distTable, distTable + supportSize + 2, 0);
      distTable[1] = distTable[1 - (supportSize + 2)]; // Get second last PMF's max.
      distTable[0] = distTable[1]; // MDR <- max
      distTable[2 + supportSize - 1] = 1; // Last prob.  
    }
  }
  
  
  // distTable[] is a space of (2 + supportSize) * Ndistr
  struct Node { double mdr, max, p0, p1, cv, *P;  };
  std::vector<Node> nodevec;
  void makeSequence(double *distTable, int Ndist, int supportSize)
    // Here, Ndist is 19001.
  {
    nodevec.resize(Ndist);
    int64_t nrow = supportSize + 2;
    cp.parFor(0, Ndist, [&](std::size_t i, std::size_t t)->bool
    {
      nodevec[i].mdr = distTable[nrow * i];
      nodevec[i].max = distTable[nrow * i + 1];
      nodevec[i].P = distTable + nrow * i + 2;
      nodevec[i].p0 = nodevec[i].P[0];
      nodevec[i].p1 = nodevec[i].P[supportSize - 1];
      if (i == 0) nodevec[i].cv = 1e300;
      else if (Ndist - 1 - i == 0) nodevec[i].cv = 0;
      else
      { 
        double delta = nodevec[i].max / (supportSize - 1);
        double ex2 = 0;
        for (int k = 0; k < supportSize; ++k)
        { 
          double a = delta * k;
          ex2 += a * a * nodevec[i].P[k];
        }
        auto &mdr = nodevec[i].mdr;
        nodevec[i].cv = std::sqrt(ex2 - mdr * mdr) / mdr;  
      } 
      return false;
    }, Ndist / (cp.maxCore * cp.maxCore * cp.maxCore) + 1); 
  }
  
  
  std::vector<int> indicesWanted;
  std::vector<Node> subseq; // Do not think this is the longest subsequence.
  
  
  int *longestMonoSeq, *longestMonoSeqEnd;
  
  
  struct P0mono {
    bool operator()(const Node &x, const Node &y) { return x.p0 > y.p0; } };
  P0mono p0mono;
  
  
  struct P0semiMono {
    bool operator()(const Node &x, const Node &y) { return x.p0 >= y.p0; } };
  P0semiMono p0semimono;
  
  
  struct P1mono {
    bool operator()(const Node &x, const Node &y) { return x.p1 < y.p1; } };
  P1mono p1mono;
  
  
  struct P1semiMono {
    bool operator()(const Node &x, const Node &y) { return x.p1 <= y.p1; } };
  P1semiMono p1semimono;
  
  
  struct CVmono {
    bool operator()(const Node &x, const Node &y) { return x.cv > y.cv; } };
  CVmono cvmono;
  
  
  struct CVsemiMono {
    bool operator()(const Node &x, const Node &y) { return x.cv >= y.cv; } };
  CVsemiMono cvsemimono;
    
  
  struct P0P1CVmono {
    bool operator()(const Node &x, const Node &y) {
      return x.p0 >= y.p0 and x.p1 <= y.p1 and x.cv >= y.cv and 
      (x.p0 != y.p0 or x.p1 != y.p1 or x.cv != y.cv); } };
  P0P1CVmono p0p1cvmono;
  
  
  struct P0P1CVsemiMono {
    bool operator()(const Node &x, const Node &y) {
      return x.p0 >= y.p0 and x.p1 <= y.p1 and x.cv >= y.cv; } };
  P0P1CVsemiMono p0p1cvsemimono;
  
  
  template <typename Lessthan>
  void longestSubseqQA(Lessthan &lthn)
    // Return a range within indicesWanted.
  {
    subseq.resize(nodevec.size());
    indicesWanted.resize(nodevec.size() * 2);
    int *C = indicesWanted.data();
    subseq[0] = nodevec[0];
    subseq.resize(1);
    auto &x = nodevec;
    for (int i = 1, iend = x.size(); i < iend; ++i)
    {
      auto it = std::lower_bound(subseq.begin(), subseq.end(), x[i], lthn);
      if (it == subseq.end())
      {
        subseq.push_back(x[i]);
        C[i] = subseq.size();
      }
      else
      {
        *it = x[i];
        C[i] = it - subseq.begin() + 1;
      }
    }
    int whichMax = std::max_element(C, C + x.size()) - C;
    int *v = &indicesWanted.back();
    int *vend = v + 1;
    *v = whichMax;
    --v;
    for (int i = whichMax - 1; i >= 0; --i)
    {
      if (C[i] == C[whichMax] - 1)
      {
        *v = i;
        whichMax = i;
        --v;
      }
    }
    longestMonoSeq = v;
    longestMonoSeqEnd = vend;
    
    
    // DO NOT DELETE! FOR CHECKING IF THE SUBSEQUENCE length.
    if constexpr (false) // 
    {
      Rcout << "longestSubseqLength() = " << longestSubseqLength() << ", ";
      Rcout << "vend - v = " << vend - v << "\n";
    }
  }
  
  
  // This function is just for checking.
  int longestSubseqLength()
  {
    auto lessthan = [](const Node &x, const Node &y)->bool
    {
      return (x.max <= y.max and x.p0 >= y.p0 and x.p1 <= y.p1 and x.cv >= y.cv)
      // and (x.max != y.max or x.p0 != y.p0 or x.p1 != y.p1 or x.cv != y.cv)
      ;
    }; 
    std::vector<int> sizes(nodevec.size(), 1);
    for (int i = 1, iend = sizes.size(); i < iend; ++i)
    {
      for (int j = i - 1; j >= 0; --j)
      {
        if ( lessthan(nodevec[j], nodevec[i]) )
          sizes[i] = std::max(sizes[i], sizes[j] + 1);
      }
    }
    return *std::max_element(sizes.begin(), sizes.end());
  }
  
  
  // Write the result in-place.
  void tuneMaxP0p1forMonotonicity() // rst[] has the same dimension as distMat.
  {
    int supportSize = nodevec[1].P - nodevec[0].P - 2;
    int llsize = longestMonoSeqEnd - longestMonoSeq;
    auto f = [&](std::size_t i, std::size_t t)->bool
    {
      int leftInd = longestMonoSeq[i - 1];
      int rightInd = longestMonoSeq[i];
      auto &right = nodevec[rightInd];
      for(int j = leftInd + 1; j < rightInd; ++j, ++leftInd)
      {
        auto &left = nodevec[leftInd];
        auto &mid = nodevec[j];
        if (leftInd == 0)
        {
          mid.max = mid.mdr / right.mdr * right.max;
          mid.P[-1] = mid.max;
          for (int u = 0; u < supportSize; ++u) mid.P[u] = right.P[u];
        }
        else
        {
          
          double leftMax = std::min(mid.mdr / left.mdr * left.max, right.max);
          double rightMax = std::max(mid.mdr / right.mdr * right.max, left.max);
          double max = (leftMax + rightMax) / 2;
          
          
          double leftMDR = max / left.max * left.mdr;
          double rightMDR = max / right.max * right.mdr;
          double w = 0;
          if (rightMDR != leftMDR) w = (mid.mdr - leftMDR) / (rightMDR - leftMDR);
          
          
          mid.max = max;
          mid.P[-1] = max; // Assign the max to the value in the actual distribution table.
          for (int u = 0; u < supportSize; ++u)
            mid.P[u] = (1 - w) * left.P[u] + w * right.P[u];
        }
      }
      return false;
    };
    
    
    auto &m = cp.maxCore;
    cp.parFor(1, llsize, f, llsize / (m * m) + 1);
  }
  
  
};



  
// TrBtable should not contain any NA, NaN or Inf.
template <int RIBlib, bool downScaleSupport>
List FTDFQAcpp(
    NumericMatrix TrBtable, 
    int *supportSizes,
    int NsupportSizes,
    String regridMethod,
    bool outputProbsInRows,
    int fineDiscretizationSize,
    int maxCore, bool verbose)
{
  // mdr, max, p0, a, b, c, d;
  int Ndist = TrBtable.ncol();
  Ftdfqa D(maxCore);
  
  
  vec<double*> ptrbank(NsupportSizes * 3);
  double **distablesResv = &ptrbank[0], **distables = &ptrbank[NsupportSizes];
  double **untuneDistables = &ptrbank[NsupportSizes * 2];
  
  
  List tunedList(NsupportSizes), untunedList(NsupportSizes);
  for (int j = 0; j < NsupportSizes; ++j)
  {
    // distablesSizes[j] = (supportSizes[j] + 2LL) * (Ndist + 2);
    NumericMatrix distTable(supportSizes[j] + 2, Ndist + 2);
    distables[j] = &distTable[0];
    distablesResv[j] = distables[j];
    tunedList[j] = distTable;
    NumericMatrix distTableBeforeTune(distTable.nrow(), distTable.ncol());
    untuneDistables[j] = &distTableBeforeTune[0];
    untunedList[j] = distTableBeforeTune;
  }
  
  
  if constexpr (false)
  {
    Rcout << "1.1\n";
    NumericMatrix dj = tunedList[0];
    Rcout << "&dj[0] = " << uint64_t(&dj[0]) << ", ";
    Rcout << "distables[0] = " << uint64_t(distables[0]) << "\n\n";
  }
  
  
  StringVector distlistNames(NsupportSizes);
  for (int j = 0; j < NsupportSizes; ++j)
    distlistNames[j] = "p" + std::to_string(supportSizes[j]);
  tunedList.attr("names") = distlistNames;
  untunedList.attr("names") = distlistNames;
  
  
  if constexpr (false)
  {
    Rcout << "1.2\n";
    NumericMatrix dj = tunedList[0];
    Rcout << "&dj[0] = " << uint64_t(&dj[0]) << ", ";
    Rcout << "distables[0] = " << uint64_t(distables[0]) << "\n\n";
  }
  
  
  if (verbose)
    Rcout << "Fine discretize TrB over " + 
      std::to_string(fineDiscretizationSize) + "-point supports...\n";
  
  
  if (regridMethod == "lr") D.discretize<0, RIBlib, downScaleSupport> (
    &TrBtable[0], Ndist, fineDiscretizationSize, &distables[0], supportSizes, NsupportSizes);
  else if (regridMethod == "lmm") D.discretize<1, RIBlib, downScaleSupport> (
    &TrBtable[0], Ndist, fineDiscretizationSize, &distables[0], supportSizes, NsupportSizes);
  else if (regridMethod == "r4") D.discretize<2, RIBlib, downScaleSupport> (
    &TrBtable[0], Ndist, fineDiscretizationSize, &distables[0], supportSizes, NsupportSizes);
  else stop("Regrid method not implemented.");
  
  
  // This is because distables's contents have been changed in D.discretize.
  std::copy(distablesResv, distablesResv + NsupportSizes, distables);
  
  
  for (int j = 0; j < NsupportSizes; ++j)
  {
    std::copy(distables[j], distables[j] + 
      (supportSizes[j] + 2LL) * (Ndist + 2LL), untuneDistables[j]);
  }
  
  
  for (int j = 0; j < NsupportSizes; ++j)
  {
    auto distTable = distables[j];
    auto supportSize = supportSizes[j];
    
    
    if (verbose) Rcout << 
      "Tune PMF table with " << supportSize << "-point supports:\n";
    
   
    D.makeSequence(&distTable[0], Ndist + 2, supportSize);
    D.longestSubseqQA( D.p0semimono ); // D.longestSubseqQA( D.p0mono );
    if (verbose)
    {
      Rcout << "Length of the longest decreasing subsequence of P0 = " << 
        D.longestMonoSeqEnd - D.longestMonoSeq << ". Tuning..\n";
    }
    D.tuneMaxP0p1forMonotonicity();
    D.makeSequence(&distTable[0], Ndist + 2, supportSize);
    
    
    // ===========================================================================
    // DO NOT DELETE! The following is supposed to print 19001 for all that have 
    //   been fixed.
    // ===========================================================================
    if constexpr (false)
    {
      D.longestSubseqQA( D.p0semimono );
      Rcout << "Size of the longest decreasing subsequence of P0 = " << 
        D.longestMonoSeqEnd - D.longestMonoSeq << "\n";
    }
    
    
    D.longestSubseqQA( D.p1semimono ); // D.longestSubseqQA( D.p1mono );
    if (verbose)
    {
      Rcout << "Length of the longest increasing subsequence of Pmax = " << 
        D.longestMonoSeqEnd - D.longestMonoSeq << ". Tuning..\n";
    }
    D.tuneMaxP0p1forMonotonicity();
    D.makeSequence(&distTable[0], Ndist + 2, supportSize);
    
    
    // ===========================================================================
    // DO NOT DELETE! The following is supposed to print 19001 for all that have 
    //   been fixed.
    // ===========================================================================
    if constexpr (false)
    {
      D.longestSubseqQA( D.p0semimono );
      Rcout << "Length of the longest decreasing subsequence of P0 = " << 
        D.longestMonoSeqEnd - D.longestMonoSeq << "\n";
      D.longestSubseqQA( D.p1semimono );
      Rcout << "Length of the longest decreasing subsequence of Pmax = " << 
        D.longestMonoSeqEnd - D.longestMonoSeq << "\n";
    }
    
    
    D.longestSubseqQA( D.cvsemimono ); // D.longestSubseqQA( D.cvmono );
    if (verbose)
    {
      Rcout << "Length of the longest decreasing subsequence of coefficients of variation = " << 
        D.longestMonoSeqEnd - D.longestMonoSeq << ". Tuning..\n";
    }
    D.tuneMaxP0p1forMonotonicity();
    D.makeSequence(&distTable[0], Ndist + 2, supportSize);
    
    
    // ===========================================================================
    // DO NOT DELETE! The following is supposed to print 19001 for all that have 
    //   been fixed.
    // ===========================================================================
    if constexpr (false)
    {
      D.longestSubseqQA( D.p0semimono );
      Rcout << "  Length of the longest decreasing subsequence of P0 = " << 
        D.longestMonoSeqEnd - D.longestMonoSeq << "\n";
      D.longestSubseqQA( D.p1semimono );
      Rcout << "  Length of the longest decreasing subsequence of Pmax = " << 
        D.longestMonoSeqEnd - D.longestMonoSeq << "\n";
      D.longestSubseqQA( D.cvsemimono );
      Rcout << "  Length of the longest decreasing subsequence of CV = " << 
        D.longestMonoSeqEnd - D.longestMonoSeq << "\n";
    }
    
    
    D.longestSubseqQA( D.p0p1cvsemimono ); // D.longestSubseqQA( D.p0p1cvmono );
    if (verbose)
    {
      Rcout << "Length of the longest fully QA compliant subsequence = " << 
        D.longestMonoSeqEnd - D.longestMonoSeq << ". Tuning..\n";
    }
    D.tuneMaxP0p1forMonotonicity();
    
    
    // Check everything.
    double maxMeanDiff = -1e300;
    double maxCVdiff = -1e300;
    double priorCv = 1e300;
    double minMaxDiff = 1e300;
    double minP1diff = 1e300;
    double maxP0diff = -1e300;
    double maxPsumDiffFrom1 = -1e300;
    double minP = 0;
    for (int64_t i = 1; i < Ndist; ++i)
    {
      double *d = &distTable[0] + i * (supportSize + 2);
      double delta = d[1] / (supportSize - 1);
      double mean = 0, ex2 = 0;
      double *p = d + 2;
      double psum = 0;
      for (int k = 0; k < supportSize; ++k)
      {
        minP = std::min(p[k], minP);
        psum += p[k];
        double kdelta = k * delta;
        double partialMean = p[k] * kdelta;
        mean += partialMean;
        ex2 += partialMean * kdelta;
      }
      maxPsumDiffFrom1 = std::max(std::abs(psum - 1), maxPsumDiffFrom1);
      maxMeanDiff = std::max(maxMeanDiff, std::abs(mean - d[0]));
      minMaxDiff = std::min(d[1] - d[1 - (supportSize + 2)], minMaxDiff);
      double sd = std::sqrt(ex2 - mean * mean);
      double cv = sd != 0 ? sd / mean : 0;
      maxCVdiff = std::max(maxCVdiff, cv - priorCv);
      priorCv = cv;
    }
    
    
    if (verbose)
    {
      auto maxele = *std::max_element( 
        D.maxIterMaxPointCntr.begin(), D.maxIterMaxPointCntr.begin() + D.cp.maxCore);
      auto maxpr = *std::max_element( 
        D.maxIterMaxPointCntr.begin() + D.cp.maxCore, D.maxIterMaxPointCntr.end());
      
      
      Rcout << "Imposing MDR took at most " << maxele << " iterations.\n";
      Rcout << "Maximum number of points reduced during tuning means = " << maxpr << "\n";  
      
      
      Rcout << "Final check:\n";
      
      
      if (maxMeanDiff < 1e-10) Rcout << "MDR == mean: True.\n";
      else Rcout << "Maximum absolute difference between MDR and mean = " << maxMeanDiff << "\n";
      if (minMaxDiff > -1e-10) Rcout << "Max monotonicity: True.\n";
      else Rcout << "Minimum difference between neighboring maxes = " << minMaxDiff << "\n";
      if (maxP0diff < 1e-10) Rcout << "P0 monotonicity: True.\n";
      else Rcout << "Maximum difference between neighboring P0s = " << maxP0diff << "\n";
      if (minP1diff > -1e-10) Rcout << "P1 monotonicity: True.\n";
      else Rcout << "Minimum difference between neighboring P1s = " << minP1diff << "\n";
      if (minP >= 0) Rcout << "Probability nonnegativity: True.\n";
      else Rcout << "Minimum probability = " << minP << "\n";
      if (maxPsumDiffFrom1 < 1e-10) Rcout << "Sum of probabilities == 1: True.\n";
      else Rcout << "Maximum absolute difference between sum of probabilities and 1 = " << 
        maxPsumDiffFrom1 << "\n";
      if (maxMeanDiff < 1e-10 and minMaxDiff > -1e-10 and maxP0diff < 1e-10 and 
            minP1diff > -1e-10 and minP >= 0 and maxPsumDiffFrom1 < 1e-10)
        Rcout << "All QA requirements have been met.\n\n";
      else warning("QA requirements have not been met.\n\n");
    }
    
    
  }
  
  
  // List tunedList(NsupportSizes), untunedList(NsupportSizes);
  if (outputProbsInRows) 
  {
    for (int j = 0; j < NsupportSizes; ++j)
    {
      NumericMatrix tmp = tunedList[j], tmp2 = untunedList[j];
      tunedList[j] = transpose(tmp);
      untunedList[j] = transpose(tmp2);
    }
  }
  return List::create(
    Named("tuned") = tunedList, 
    Named("untuned") = untunedList);
}




// 7/23 night stopped here.
template <int RIBlib>
List FTDFQAcpp(
    NumericMatrix TrBtable, 
    int* supportSizes, int NsupportSizes,
    String regridMethod,
    bool outputProbsInRows,
    int fineDiscretizationSize,
    int maxCore, bool verbose, bool downScaleSupport)
{
  if (downScaleSupport)
    return FTDFQAcpp<RIBlib, true>(
      TrBtable, supportSizes, NsupportSizes,
      regridMethod,
      outputProbsInRows,
      fineDiscretizationSize,
      maxCore, verbose);
  return FTDFQAcpp<RIBlib, false>(
      TrBtable, supportSizes, NsupportSizes,
      regridMethod,
      outputProbsInRows,
      fineDiscretizationSize,
      maxCore, verbose);
}






//'
//' Inflated Transformed Beta Discretization
//' 
//' \strong{F}ine-tuned \strong{D}iscretization for \strong{F}ulfilling \strong{QA}
//' requirements, i.e. FTDFQA. The function produces a PMF table from an 
//' inflated TrB parameter table.
//' 
//' @param TrBtable  A numeric matrix of 7 rows. Each column is a vector of
//' \code{(MDR, max, }\eqn{P_0}\code{, a, b, c, d)}.
//' 
//' @param supportSizes  An integer vector, e.g. \code{c(64, 42)}.
//' 
//' @param regridMethod  Regrid method: \code{"lr"} (linear regrid) or 
//' \code{"lmm"} (local moment matching) or \code{"r4"} (four-point regrid).
//' Default to \code{"lr"}.
//' 
//' @param outputProbsInRows  A boolean. \code{TRUE} will produce a PMF 
//' table (matrix) where every row is a vector of (\code{MDR, max}, \eqn{P_0}, 
//' \code{P1, ..., Pmax}). \code{FALSE} will store it in a column. 
//' Default to \code{FALSE}.
//' 
//' @param fineDiscretizationSize  An integer. Support size for raw 
//' discretization before tunning. Default to 2000.
//' 
//' @param verbose  A boolean. \code{TRUE} prints QA checking messages.
//' 
//' @param downScaleSupport  A boolean. When the raw discretization's mean is
//' less than MDR, \code{TRUE} will downscale the support to match the MDR, 
//' while \code{FALSE} will increase \eqn{P_0} and downscale the main part 
//' probabilities to match the MDR. Default to \code{FALSE}, which is the 
//' better choice from experiments.
//' 
//' @inheritParams LBFGSBtrbFitList
//' 
//' @inherit  computeWindows details
//' 
//' @return A list of two named objects:
//' \itemize{
//' \item\code{distTable}: A list of final PMF tables corresponding to 
//' \code{supportSizes}. 
//' 
//' Names of the tables in \code{distTable} are \code{p} succeeded by 
//' the support sizes. For example, 
//' if \code{supportSizes = c(42L, 64L)}, then the names would be
//' \code{'p42'} and \code{'p64'}.
//' 
//' If \code{outputProbsInRows = FALSE}, every column 
//' of \code{distTable$p42} or \code{distTable$p64} would be 
//' \code{MDR, max, P0, P1, ..., Pmax}. Otherwise the tables are transposed.
//' 
//' \item\code{distTableBeforeTune}: A list of PMF tables before tuning for
//' monotonicity requirements. These are returned for diagnostics.
//' }
//' 
//' @example inst/examples/FTDFQA.R
//' 
//' 
//' 
//' 
// [[Rcpp::export]]
List FTDFQA(
    NumericMatrix TrBtable, 
    IntegerVector supportSizes, 
    String regridMethod = "lr",
    String RIBlib = "Numerical Recipes",
    bool outputProbsInRows = false,
    int fineDiscretizationSize = 2000,
    int maxCore = 1000, 
    bool verbose = true, 
    bool downScaleSupport = false)
{
  fineDiscretizationSize = std::max(
    *std::max_element(supportSizes.begin(), supportSizes.end()) * 2, 
    fineDiscretizationSize);
  int NsupportSizes = supportSizes.size();
  
  
  if (RIBlib == "Boost") 
    return FTDFQAcpp<0> (TrBtable, &supportSizes[0], NsupportSizes, regridMethod, outputProbsInRows,
                         fineDiscretizationSize, maxCore, verbose, downScaleSupport);
  if (RIBlib == "R::pbeta")
    return FTDFQAcpp<1> (TrBtable, &supportSizes[0], NsupportSizes, regridMethod, outputProbsInRows,
                         fineDiscretizationSize, maxCore, verbose, downScaleSupport);
  if (RIBlib == "Numerical Recipes")
    return FTDFQAcpp<2> (TrBtable, &supportSizes[0], NsupportSizes, regridMethod, outputProbsInRows,
                         fineDiscretizationSize, maxCore, verbose, downScaleSupport);
  stop("Unknown Regularized Incomplete Beta function.");
  return List::create();
}
  

  

#undef vec








