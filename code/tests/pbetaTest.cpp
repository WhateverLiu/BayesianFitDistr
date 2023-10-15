// [[Rcpp::plugins(cpp17)]]
#include <boost/math/special_functions/beta.hpp>
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;


#include "../cpp/h/textbookRegularizedIncompleteBeta.h"


// Compute the regularized incomplete Beta function.
// [[Rcpp::export]]
NumericVector RIBF(NumericVector q, NumericVector a, NumericVector b, 
                   int choice, double fpmin, double eps)
{
  NumericVector rst(q.size());
  using namespace boost::math;
  using namespace boost::math::policies;
  for (int i = 0, iend = q.size(); i < iend; ++i)
  {
    if (choice == 0) rst[i] = boost::math::ibeta( a[i], b[i], q[i], 
        policy<digits10<10> >()
      );
    else if (choice == 1) rst[i] = R::pbeta( q[i], a[i], b[i], 1, 0 );
    else rst[i] = textbookRIB::betai(q[i], a[i], b[i], fpmin, eps);
  }
  return rst;
}
































