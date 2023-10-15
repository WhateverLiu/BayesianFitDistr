#pragma once
// #include <boost/math/special_functions/beta.hpp>


// Regularized incomplete beta function.
template<int RIBlib>
inline double ribetaf(
    double u, double a, double b, int lowerTail, int giveLog)
{
  double rst = 0;
  if (RIBlib == 0)
  {
    rst = lowerTail ? boost::math::ibeta(a, b, u) : 
      boost::math::ibetac(a, b, u);
    rst = giveLog ? std::log(rst) : rst;
  }
  else if (RIBlib == 1)
  {
    rst = R::pbeta(u, a, b, lowerTail, giveLog);
  }
  else if (RIBlib == 2)
  {
    rst = NumericalRecipesRIB::betai(u, a, b, 1e-38, 1e-13);
    rst = lowerTail ? rst : 1.0 - rst;
    rst = giveLog ? std::log(rst) : rst;
  }
  return rst;
}


template<int RIBlib>
inline double ribetaf_inv(
    double u, double a, double b, int lowerTail, int giveLog)
{
  double rst = 0;
  if (RIBlib == 0)
  { 
    rst = lowerTail ? boost::math::ibeta_inv(a, b, u) : 
      boost::math::ibetac_inv(a, b, u);
    rst = giveLog ? std::log(rst) : rst;
  } 
  else if (RIBlib == 1)
  {
    rst = R::qbeta(u, a, b, lowerTail, giveLog);
  }
  else if (RIBlib == 2)
  {
    u = lowerTail ? u : 1.0 - u;
    rst = NumericalRecipesRIB::invbetai(u, a, b, 1e-38, 1e-13); 
    rst = giveLog ? std::log(rst) : rst;
  }
  return rst;
} 




// =============================================================================
// Made for tests.
// =============================================================================
inline double ribetaf(
    double u, double a, double b, int lowerTail, int giveLog, bool useBoost)
{
  double rst = 0;
  if (useBoost)
  {
    rst = lowerTail ? boost::math::ibeta(a, b, u) : boost::math::ibetac(a, b, u);
    rst = giveLog ? std::log(rst) : rst;
  }
  else
  {
    rst = R::pbeta(u, a, b, lowerTail, giveLog);
  }
  return rst;
}


// template<bool useBoost>
inline double ribetaf_inv(
    double u, double a, double b, int lowerTail, int giveLog, bool useBoost)
{
  double rst = 0;
  if (useBoost)
  { 
    rst = lowerTail ? boost::math::ibeta_inv(a, b, u) : 
    boost::math::ibetac_inv(a, b, u);
    rst = giveLog ? std::log(rst) : rst;
  } 
  else
  {
    rst = R::qbeta(u, a, b, lowerTail, giveLog);
  }
  return rst;
} 































