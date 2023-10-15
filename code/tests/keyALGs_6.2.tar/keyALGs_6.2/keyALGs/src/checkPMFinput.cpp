// #pragma once
#include <Rcpp.h>


// Check ordering and duplicate of input PMF's support points
void checkPMFintegrity(Rcpp::List X, bool checkReverseOrder = false)
{
  Rcpp::NumericVector val = X[0], P = X[1];
  if (P[0] < 0) Rcpp::warning("Input PMF has negative probability. Program could crash. Result can be wrong.");
  double psum = P[0];
  bool earlyBreak = false;
  for (int i = 1, iend = val.size(); i < iend; ++i)
  {
    bool increasing = val[i - 1] < val[i];
    bool decreasing = val[i - 1] > val[i];
    if (  !(  (increasing and !checkReverseOrder) or (decreasing and checkReverseOrder)  )  )
    {
      Rcpp::warning("Input PMF's support is unsorted or has duplicate points. Program could crash. Result can be wrong.");
      earlyBreak = true;
      break;
    }
    if (P[i] < 0)
    {
      Rcpp::warning("Input PMF has negative probability. Program could crash. Result can be wrong.");
      earlyBreak = true;
      break;
    }
    psum += P[i];
  }
  if (!earlyBreak and std::abs(psum - 1.0) > 1e-5)
    Rcpp::warning("Input PMF probability sum < 1 - 1e-5 OR > 1 + 1e-5.");
}


void checkPMFlistIntegrity(Rcpp::List distlist, bool checkReverseOrder = false)
{
  for (int j = 0, jend = distlist.size(); j < jend; ++j)
  {
    Rcpp::List X = distlist[j];
    checkPMFintegrity(X, checkReverseOrder);
  }
}



































