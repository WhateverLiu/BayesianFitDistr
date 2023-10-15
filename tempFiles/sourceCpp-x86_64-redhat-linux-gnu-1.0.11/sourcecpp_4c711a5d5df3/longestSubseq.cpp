// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#define vec std::vector


// Return the one-based indices.
// [[Rcpp::export]]
IntegerVector longestIncreasingSubseqSlow(NumericVector x)
{
  std::vector<int> v(x.size() * 2, 0);
  int *C = &v[0];
  C[0] = 1;
  // Find the size of the longest increasing subsequence ending at x[i].
  // The size is stored in C[i].
  for (int i = 1, iend = x.size(); i < iend; ++i)
  {
    for (int j = i - 1; j >= 0; --j)
    {
      if (x[i] > x[j])
        C[i] = std::max(C[i], C[j] + 1);
    }
  }
  v.resize(x.size());
  int whichMax = std::max_element(C, C + x.size()) - C;
  v.push_back(whichMax + 1);
  for (int i = whichMax - 1; i >= 0; --i)
  {
    if (C[i] == C[whichMax] - 1)
    {
      v.push_back(i + 1);
      whichMax = i;
    }
  }
  std::reverse(v.begin() + x.size(), v.end());
  return IntegerVector(v.begin() + x.size(), v.end());
}




// [[Rcpp::export]]
IntegerVector longestIncreasingSubseq(NumericVector x) 
{
  std::vector<int> v(x.size() * 2, 0);
  int *C = &v[0];
  C[0] = 1;
  // Find the size of the longest increasing subsequence ending at x[i].
  // The size is stored in C[i].
  std::vector<double> subseq(x.size());
  subseq.resize(1);
  subseq[0] = x[0];
  for (int i = 1, iend = x.size(); i < iend; ++i)
  {
    auto it = std::lower_bound(subseq.begin(), subseq.end(), x[i]);
    if (it >= subseq.end())
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
  v.resize(x.size());
  int whichMax = std::max_element(C, C + x.size()) - C;
  v.push_back(whichMax + 1);
  for (int i = whichMax - 1; i >= 0; --i)
  { 
    if (C[i] == C[whichMax] - 1)
    { 
      v.push_back(i + 1);
      whichMax = i;
    } 
  } 
  std::reverse(v.begin() + x.size(), v.end());
  return IntegerVector(v.begin() + x.size(), v.end());
} 




//'  
//' Longest Nondecreasing Subsequence.
//' 
//' Compute the longest nondecreasing subsequence from a given sequence of numbers.
//' 
//' @param x  A numeric vector.
//' 
//' @return An integer vector of the 1-based indices of the elements in the 
//' subsequence of interest.
//' 
//' @details More details can be found \href{https://en.wikipedia.org/wiki/Longest_increasing_subsequence}{here}.
//' 
//' @example inst/examples/longestNonDecreasingSubseq.R
//'
// [[Rcpp::export]]
IntegerVector longestNonDecreasingSubseq(NumericVector x) 
{
  std::vector<int> v(x.size() * 2, 0);
  int *C = &v[0];
  C[0] = 1;
  // Find the size of the longest increasing subsequence ending at x[i].
  // The size is stored in C[i].
  std::vector<double> subseq(x.size());
  subseq.resize(1);
  subseq[0] = x[0];
  for (int i = 1, iend = x.size(); i < iend; ++i)
  {
    // auto it = std::upper_bound(subseq.begin(), subseq.end(), x[i]);
    auto it = std::lower_bound(
      subseq.begin(), subseq.end(), x[i], [](const double &u, const double &v)->bool
        {
          return u <= v;
        });
    if (it >= subseq.end())
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
  v.resize(x.size());
  int whichMax = std::max_element(C, C + x.size()) - C;
  v.push_back(whichMax + 1);
  for (int i = whichMax - 1; i >= 0; --i)
  { 
    if (C[i] == C[whichMax] - 1)
    { 
      v.push_back(i + 1);
      whichMax = i;
    } 
  } 
  std::reverse(v.begin() + x.size(), v.end());
  return IntegerVector(v.begin() + x.size(), v.end());
} 




#undef vec
















#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// longestIncreasingSubseqSlow
IntegerVector longestIncreasingSubseqSlow(NumericVector x);
RcppExport SEXP sourceCpp_56_longestIncreasingSubseqSlow(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(longestIncreasingSubseqSlow(x));
    return rcpp_result_gen;
END_RCPP
}
// longestIncreasingSubseq
IntegerVector longestIncreasingSubseq(NumericVector x);
RcppExport SEXP sourceCpp_56_longestIncreasingSubseq(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(longestIncreasingSubseq(x));
    return rcpp_result_gen;
END_RCPP
}
// longestNonDecreasingSubseq
IntegerVector longestNonDecreasingSubseq(NumericVector x);
RcppExport SEXP sourceCpp_56_longestNonDecreasingSubseq(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(longestNonDecreasingSubseq(x));
    return rcpp_result_gen;
END_RCPP
}
