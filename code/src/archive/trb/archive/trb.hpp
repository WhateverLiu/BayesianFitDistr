// #include <R.h>
// #include <Rmath.h>


// #define ACT_DLIM__0(x, y)   (R_FINITE(x) ? R_pow(x, y) : 0.)


double betaint_raw(double x, double a, double b, double x1m)
{
  /* Here, assume that (x, a, b) are not NA, 0 < x < 1 and 0 < a < Inf. */
  
  if (b > 0)
  {
    /* I(x, a, b) = 1 - I(1 - x, b, a) */
    double Ix = (x > 0.5) ?
    // pbeta(x1m, b, a, /*l._t.*/0, /*give_log*/0) :
    R::pbeta(x1m, b, a, /*l._t.*/0, /*give_log*/0) :
    // pbeta(x,   a, b, /*l._t.*/1, /*give_log*/0);
    R::pbeta(x,   a, b, /*l._t.*/1, /*give_log*/0);
    // return gammafn(a) * gammafn(b) * Ix;
    return std::tgamma(a) * std::tgamma(b) * Ix;
  }
  
  // double r = floor(-b);
  double r = std::floor(-b);
  
  // ACT_nonint(x):  fabs( x - ACT_forceint(x) ) > 1e-7 * fmax2(1., fabs(x) )
  auto notInt = [](double x)->bool
  {
    return std::abs( x - int(x) ) > 1e-7 * std::max(1.0, std::abs(x) );
  };
  
  
  // if (! (ACT_nonint(b) && a - r - 1 > 0)) return R_NaN;
  if ( ! ( notInt(b) and a - r - 1 > 0 ) ) 
  {
    // Rcpp::Rcout << "!!\n";
    return 1e300; 
  }
  
  
  /* There are two quantities to accumulate in order to compute the
   * final result: the alternating sum (to be stored in 'sum') and
   * the ratio [(a - 1) ... (a - r)]/[b(b + 1) ... (b + r)] (to be
   * stored in 'ratio'). Some calculations are done in the log
   * scale. */
  int i;
  double ap = a, bp = b;	/* copies of a and b */
  // double lx = log(x);		/* log(x) */
  double lx = std::log(x);		/* log(x) */
  // double lx1m = log(x1m);	/* log(1 - x) */
  double lx1m = std::log(x1m);	/* log(1 - x) */
  // double x1 = exp(lx1m - lx);	/* (1 - x)/x */
  double x1 = std::exp(lx1m - lx);	/* (1 - x)/x */
  double c, tmp, sum, ratio;
  
  /* Computation of the first term in the alternating sum. */
  ap--;			     /* a - 1 */
  // c = exp(ap * lx + bp * lx1m)/bp; /* (x^(a - 1) (1 - x)^b) / b */
  c = std::exp(ap * lx + bp * lx1m)/bp; /* (x^(a - 1) (1 - x)^b) / b */
  sum = c;			     /* first term */
  ratio = 1/bp;		     /* 1 / b */
  bp++;			     /* b + 1 */
  
  /* Other terms in the alternating sum iff r > 0.
   * Relies on the fact that each new term in the sum is
   *
   *  [previous term] * (a - i - 1)(1 - x)/[(b + i + 1) x]
   *
   * for i = 0, ..., r - 1. We need to compute this value as
   *
   *  {[previous term] * [(1 - x)/x]} * [(a - i - 1)/(b + i + 1)]
   *
   * to preserve accuracy for very small values of x (near
   * DBL_MIN).
   */
  for (i = 0; i < r; i++)
  {
    tmp = ap/bp;		/* (a - i - 1)/(b + i + 1) */
  c = tmp * (c * x1);	/* new term in the sum  */
  sum += c;
  ratio *= tmp;
  ap--;
  bp++;
  }
  
  /* I(x, a, b) = 1 - I(1 - x, b, a) */
  double lIx = (x > 0.5) ?
  R::pbeta(x1m, bp, ap, /*l._t.*/0, /*give_log*/1) :
    R::pbeta(x,   ap, bp, /*l._t.*/1, /*give_log*/1);
  
  // return(-gammafn(a + b) * sum + (ratio * ap) * exp(lgammafn(ap) + lgammafn(bp) + lIx));
  
  return -std::tgamma(a + b) * sum + (ratio * ap) * 
    exp(std::lgamma(ap) + std::lgamma(bp) + lIx);
  
}


double levtrbeta(// double limit, 
                 double shape1, double shape2, double shape3,
                 double scale //, double order, int give_log
                   )
{
// #ifdef IEEE_754
//   if (ISNAN(limit) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) ||
//       ISNAN(scale) || ISNAN(order))
//     return limit + shape1 + shape2 + shape3 + scale + order;
// #endif
//   if (!R_FINITE(shape1) ||
//       !R_FINITE(shape2) ||
//       !R_FINITE(shape3) ||
//       !R_FINITE(scale) ||
//       !R_FINITE(order) ||
//       shape1 <= 0.0 ||
//       shape2 <= 0.0 ||
//       shape3 <= 0.0 ||
//       scale  <= 0.0)
//     return R_NaN;
  
  // if (order <= - shape3 * shape2) return R_PosInf;
  
  // if (limit <= 0.0) return 0.0;
  // constexpr const double order = 1;
  
  
  double logv, u, u1m, Ix;
  // double tmp = order / shape2;
  double tmp = 1.0 / shape2;
  
  
  // logv = shape2 * (std::log(limit) - std::log(scale));
  logv = shape2 * ( - std::log(scale));
  u = std::exp( -log1pexp(-logv) );
  u1m = std::exp( -log1pexp(logv) );
  
  
  Ix = (u > 0.5) ?
  R::pbeta(u1m, shape1, shape3, /*l._t.*/1, /*give_log*/0) :
    R::pbeta(u,   shape3, shape1, /*l._t.*/0, /*give_log*/0);
  
  
  return 
    // R_pow(scale, order)
    scale
    * betaint_raw(u, shape3 + tmp, shape1 - tmp, u1m)
    // / (gammafn(shape1) * gammafn(shape3))
    / (std::tgamma(shape1) * std::tgamma(shape3))
    // + ACT_DLIM__0(limit, order) * Ix;
    + Ix;

}
// #define ACT_DLIM__0(x, y)   (R_FINITE(x) ? R_pow(x, y) : 0.)


double dtrbeta(double x, double shape1, double shape2, double shape3,
               double scale, bool give_log)
{
// #ifdef IEEE_754
//   if (ISNAN(x) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) || ISNAN(scale))
//     return x + shape1 + shape2 + shape3 + scale;
// #endif
//   if (!R_FINITE(shape1) ||
//       !R_FINITE(shape2) ||
//       !R_FINITE(shape3) ||
//       shape1 <= 0.0 ||
//       shape2 <= 0.0 ||
//       shape3 <= 0.0 ||
//       scale <= 0.0)
//     return R_NaN;

  // if (!R_FINITE(x) || x < 0.0)
  //   return ACT_D__0;


  /* handle x == 0 separately */
  if (x == 0.0)
  {
    // if (shape2 * shape3 < 1) return R_PosInf;
    if (shape2 * shape3 < 1) return 1e300;
    // if (shape2 * shape3 > 1) return ACT_D__0;
    if (shape2 * shape3 > 1) return 0;
    /* else */
    // return give_log ?
    // log(shape2) - log(scale) - lbeta(shape3, shape1) :
    //   shape2 / (scale * beta(shape3, shape1));
    return give_log ?
    std::log(shape2) - std::log(scale) - R::lbeta(shape3, shape1) :
      shape2 / (scale * std::beta(shape3, shape1));
  }

  double logv, logu, log1mu;

  // logv = shape2 * (log(x) - log(scale));
  logv = shape2 * (std::log(x) - std::log(scale));
  logu = - log1pexp(-logv);
  log1mu = - log1pexp(logv);

  // return ACT_D_exp( log(shape2) + shape3 * logu + shape1 * log1mu
  //                    - log(x) - lbeta(shape3, shape1) );
  double logP = std::log(shape2) + shape3 * logu + shape1 * log1mu
    -std::log(x) - R::lbeta(shape3, shape1);
  return give_log ? logP : std::exp(logP);
}

// #define ACT_D_exp(x)      (log_p  ?  (x)   : exp(x))




double ptrbeta(double q, double shape1, double shape2, double shape3,
               double scale, int lower_tail, int log_p)
{
// #ifdef IEEE_754
//   if (ISNAN(q) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) || ISNAN(scale))
//     return q + shape1 + shape2 + shape3 + scale;
// #endif
//   if (!R_FINITE(shape1) ||
//       !R_FINITE(shape2) ||
//       !R_FINITE(shape3) ||
//       shape1 <= 0.0 ||
//       shape2 <= 0.0 ||
//       shape3 <= 0.0 ||
//       scale  <= 0.0)
//     return R_NaN;
//   
  if (q <= 0) return 0;
  
  double logvm, u;
  
  logvm = shape2 * (std::log(scale) - std::log(q)); /* -log v */
    u = std::exp(-log1pexp(logvm));
  
    if (u > 0.5)
    {
      /* Compute (1 - x) accurately */
      double u1m = std::exp(-log1pexp(-logvm));
      return R::pbeta(u1m, shape1, shape3, 1 - lower_tail, log_p);
    }
    
    /* else u <= 0.5 */
    return R::pbeta(u, shape3, shape1, lower_tail, log_p);
}




double trbLimMean(double a, double b, double c, double d)
{
  return levtrbeta(a / c, c, b / c, d);
}


double trbLogPdf(double x, double a, double b, double c, double d)
{
  dtrbeta(x, a / c, c, b / c, d, true);
}


double trbCdf(double x, double a, double b, double c, double d)
{
  
}















































