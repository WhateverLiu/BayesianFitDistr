// [[Rcpp::plugins(cpp11)]]
# include <Rcpp.h>
using namespace Rcpp;


// P[size - 1] != 0 and P[0] != 0 have been guaranteed.
template<bool redistribute>
double densityWhenXwithin(double *val, double *P, int size, double x)
{
  int k = std::lower_bound(val, val + size, x) - val;
  if(k == 0) return P[0] / (val[1] - val[0]);
  double rst = 0;
  int right = k, left = k - 1;
  while(P[right] == 0) ++right;
  while(P[left] == 0) --left;
  double intervalInv = 1.0 / (val[right] - val[left]);
  if(redistribute)
  {
    double w = (x - val[left]) * intervalInv;
    rst = P[left] * (1 - w) + P[right] * w;
  }
  else
  {
    if(x - val[left] <= val[right] - x) rst = P[left];
    else rst = P[right];
  }
  return rst * intervalInv;
}


template<bool redistribute>
double score(double *val, double *P, int size, double x)
{
  double ignoranceScore = 0;
  double *valend = val + size, *Pend = P + size;
  while(val < valend and P[0] == 0) { ++val; ++P; }
  while(valend > val and Pend[-1] == 0) { --valend; --Pend; }
  size = Pend - P;
  constexpr const double log2exp1 = std::log2(std::exp(1));
  
  
  if(x < val[0])
  {
    double lambda = -std::log(P[0]) / val[0]; // (x + y) / 2 = (val[0] + val[size - 1]) / 2
    // double y = val[0] + val[size - 1] - x;
    double y = val[0] * 2 - x;
    ignoranceScore = std::log(lambda) - lambda * y;
    ignoranceScore *= log2exp1;
  }
  else if(x > val[size - 1])
  {
    double lambda = -std::log(P[size - 1]) / val[size - 1];
    ignoranceScore = std::log(lambda) - lambda * x;
    ignoranceScore *= log2exp1;
  }
  else
  {
    double density = densityWhenXwithin<redistribute> (val, P, size, x);
    ignoranceScore = std::log2(density);
  }
  return -ignoranceScore;
}


// [[Rcpp::export]]
NumericVector igscore(List X, NumericVector x, bool redistributePwithin = true)
{
  NumericVector valv = X[0], Pv = X[1];
  double *val = &valv[0], *P = &Pv[0];
  int size = valv.size();
  NumericVector scores(x.size());
  if(redistributePwithin)
  {
    for(int i = 0, iend = x.size(); i < iend; ++i)
      scores[i] = score<true> (val, P, size, x[i]);
  }
  else
  {
    for(int i = 0, iend = x.size(); i < iend; ++i)
      scores[i] = score<false> (val, P, size, x[i]);
  }
  return scores;
}


// [[Rcpp::export]]
double igscoreMean(List X, List truePMF, 
                          bool redistributePwithin = true, 
                          bool symmetricScore = false)
{
  // (fx, fp) is the candidate PMF. (px, pp) is the true PMF.
  NumericVector fx = X[0], fp = X[1], px= truePMF[0], pp = truePMF[1];
  
  // double *val = &valv[0], *P = &Pv[0];
  // int size = valv.size();
  // NumericVector scores(x.size());
  double rst = 0;
  if(redistributePwithin)
  {
    for(int i = 0, iend = px.size(); i < iend; ++i)
      rst += score<true> (&fx[0], &fp[0], fx.size(), px[i]) * pp[i];
    if(symmetricScore)
    {
      for(int i = 0, iend = fx.size(); i < iend; ++i)
        rst += score<true> (&px[0], &pp[0], px.size(), fx[i]) * fp[i];
    }
  }
  else
  {
    for(int i = 0, iend = px.size(); i < iend; ++i)
      rst += score<false> (&fx[0], &fp[0], fx.size(), px[i]) * pp[i];
    if(symmetricScore)
    {
      for(int i = 0, iend = fx.size(); i < iend; ++i)
        rst += score<false> (&px[0], &pp[0], px.size(), fx[i]) * fp[i];
    }
  }
  if(symmetricScore) rst /= 2.0;
  return rst;
}











