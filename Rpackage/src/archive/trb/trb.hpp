

double betaint_raw(double x, double a, double b, double x1m)
{

  
  if (b > 0)
  {
    double Ix = (x > 0.5) ?
    R::pbeta(x1m, b, a, 0, 0) :
    R::pbeta(x,   a, b, 1, 0);
    return std::tgamma(a) * std::tgamma(b) * Ix;
  }
  

  double r = std::floor(-b);
  auto notInt = [](double x)->bool
  {
    return std::abs( x - int(x) ) > 1e-7 * std::max(1.0, std::abs(x) );
  };
  
  
  if ( ! ( notInt(b) and a - r - 1 > 0 ) ) { return 1e300; }
  
  
  int i;
  double ap = a, bp = b;
  double lx = std::log(x);
  double lx1m = std::log(x1m);
  double x1 = std::exp(lx1m - lx);
  double c, tmp, sum, ratio;
  
  
  ap--;
  c = std::exp(ap * lx + bp * lx1m) / bp;
  sum = c;
  ratio = 1/bp;
  bp++;
  
  
  for (i = 0; i < r; i++)
  {
    tmp = ap/bp;		
    c = tmp * (c * x1);	
    sum += c;
    ratio *= tmp;
    ap--;
    bp++;
  }
  
  
  double lIx = (x > 0.5) ? R::pbeta(x1m, bp, ap, 0, 1) : 
    R::pbeta(x, ap, bp, 1, 1);
  
  
  return -std::tgamma(a + b) * sum + (ratio * ap) * 
    exp(std::lgamma(ap) + std::lgamma(bp) + lIx);
}


struct BetaIntRaw
{
  double a, b;
  double tgamma_a, tgamma_b;
  double r;
  double tgamma_a_plus_b;
  double apFinal, bpFinal;
  double lgamma_ap, lgamma_bp;
  bool earlyReturn;
  
  
  BetaIntRaw(){}
  BetaIntRaw(double a, double b) { this->a = a; this->b = b; }
  
  
  bool notInt(double x)
  {
    return std::abs( x - int(x) ) > 1e-7 * std::max(1.0, std::abs(x) );
  };
  
  
  void reset(double a, double b)
  {
    this->a = a; this->b = b;
    tgamma_a = tgamma_b = 0; 
    if (b > 0)
    {
      tgamma_a = std::tgamma(a);
      tgamma_b = std::tgamma(b);
    }
    r = std::floor(-b);
    // earlyReturn = ! ( notInt(b) and a - r - 1 > 0 );
    earlyReturn = !notInt(b) or a - r - 1 <= 0 ;
    tgamma_a_plus_b = lgamma_ap = lgamma_bp = 0;
    if (earlyReturn) return;
    tgamma_a_plus_b = std::tgamma(a + b);
    apFinal = a - 1 - r; // apFinal = apFinal - 1 - r;
    bpFinal = b + 1 + r; // bpFinal = bpFinal + 1 + r;
    lgamma_ap = std::lgamma(apFinal);
    lgamma_bp = std::lgamma(bpFinal);
  }
  
  
  double operator()(double x, double x1m)
  {
    if (b > 0)
    {
      double Ix = (x > 0.5) ?
      R::pbeta(x1m, b, a, 0, 0) :
      R::pbeta(x,   a, b, 1, 0);
      return tgamma_a * tgamma_b * Ix;
    }
    
    
    // if ( ! ( notInt(b) and a - r - 1 > 0 ) ) { return 1e300; }
    if ( earlyReturn ) { return 1e300; }
    
    
    int i;
    double ap = a, bp = b;
    double lx = std::log(x);
    double lx1m = std::log(x1m);
    double x1 = std::exp(lx1m - lx);
    double c, tmp, sum, ratio;
    
    
    ap--;
    c = std::exp(ap * lx + bp * lx1m) / bp;
    sum = c;
    ratio = 1 / bp;
    bp++;
    
    
    for (i = 0; i < r; i++)
    {
      tmp = ap / bp;		
      c = tmp * (c * x1);	
      sum += c;
      ratio *= tmp;
      ap--;
      bp++;
    }
    
    
    double lIx = (x > 0.5) ? R::pbeta(x1m, bp, ap, 0, 1) : 
      R::pbeta(x, ap, bp, 1, 1);
    
    
    return -tgamma_a_plus_b * sum + (ratio * ap) * 
      exp(lgamma_ap + lgamma_bp + lIx);
  }
};


// =============================================================================
// Just for reference
// =============================================================================
double levtrbeta(double shape1, double shape2, double shape3, double scale)
{
  double logv, u, u1m, Ix;
  double tmp = 1.0 / shape2;
  logv = shape2 * ( - std::log(scale));
  u = std::exp( -log1pexp(-logv) );
  u1m = std::exp( -log1pexp(logv) );
  
  
  Rcpp::Rcout << u << ", " << u1m << "\n";
  
  
  Ix = (u > 0.5) ?
  R::pbeta(u1m, shape1, shape3, 1, 0) :
    R::pbeta(u,   shape3, shape1, 0, 0);
  
  
  return 
    scale * betaint_raw(u, shape3 + tmp, shape1 - tmp, u1m)
    / (std::tgamma(shape1) * std::tgamma(shape3)) + Ix;
}




struct TrbLimMean
{
  double shape1, shape2, shape3, scale;
  double tmp;
  double tgamma_shape1_x_tgamma_shape3;
  BetaIntRaw betaint_raw;
  
  
  void reset(double a, double b, double c)
  {
    shape1 = a / c; shape2 = c; shape3 = b / c;
    tmp = 1.0 / shape2;
    tgamma_shape1_x_tgamma_shape3 = std::tgamma(shape1) * std::tgamma(shape3);
    betaint_raw.reset(shape3 + tmp, shape1 - tmp);
  }
  
  
  double operator()(double d)
  {
    scale = d;
    double logv, u, u1m, Ix;
    logv = shape2 * ( - std::log(scale));
    u = std::exp( -log1pexp(-logv) );
    u1m = std::exp( -log1pexp(logv) );
    
    
    Ix = (u > 0.5) ?
    R::pbeta(u1m, shape1, shape3, 1, 0) :
      R::pbeta(u, shape3, shape1, 0, 0);
    
    
    // return scale * betaint_raw(
    //     u, shape3 + tmp, shape1 - tmp, u1m) / 
    //       tgamma_shape1_x_tgamma_shape3 + Ix;
    return scale * betaint_raw(u, u1m) / 
      tgamma_shape1_x_tgamma_shape3 + Ix;
  }
  
  
  
  
};











// =============================================================================
// Just for reference.
// =============================================================================
double dtrbeta(double x, double shape1, double shape2, double shape3,
               double scale, bool give_log)
{
  if (x == 0.0)
  {
    if (shape2 * shape3 < 1) return 1e300;
    if (shape2 * shape3 > 1) return 0;
    return give_log ? 
    std::log(shape2) - std::log(scale) - R::lbeta(shape3, shape1) : 
      shape2 / (scale * std::beta(shape3, shape1));
  }

  double logv, logu, log1mu;
  logv = shape2 * (std::log(x) - std::log(scale));
  logu = - log1pexp(-logv);
  log1mu = - log1pexp(logv);
  
  
  double logP = std::log(shape2) + shape3 * logu + shape1 * log1mu
    -std::log(x) - R::lbeta(shape3, shape1);
  return give_log ? logP : std::exp(logP);
}




struct TrbLogPdf
{
  double shape1, shape2, shape3, scale;
  double shape2_x_shape3;
  double log_shape2;
  double log_scale;
  double lbeta_shape3_shape1;
  double beta_shape3_shape1;
  double x0giveLog;
  double x0giveP;
  
  
  TrbLogPdf(){}
  TrbLogPdf(double a, double b, double c, double d) { reset(a, b, c, d); }
  void reset(double a, double b, double c, double d)
  {
    shape1 = a / c;
    shape2 = c;
    shape3 = b / c;
    scale  = d;
    shape2_x_shape3 = shape2 * shape3;
    log_shape2 = std::log(shape2);
    log_scale = std::log(scale);
    lbeta_shape3_shape1 = R::lbeta(shape3, shape1);
    beta_shape3_shape1 = std::beta(shape3, shape1);
    x0giveLog = log_shape2 - log_scale - lbeta_shape3_shape1;
    x0giveP = shape2 / (scale * beta_shape3_shape1);
  }
  
  
  double operator()(double x, bool give_log = true)
  {
    if (x == 0.0)
    {
      if (shape2_x_shape3 < 1) return 1e300;
      if (shape2_x_shape3 > 1) return 0;
      return give_log ? x0giveLog : x0giveP;
    }
    
    double logv, logu, log1mu;
    double logx = std::log(x);
    // logv = shape2 * (logx - std::log(scale));
    logv = shape2 * (logx - log_scale);
    logu = - log1pexp(-logv);
    log1mu = - log1pexp(logv);
    
    
    // double logP = std::log(shape2) + shape3 * logu + shape1 * log1mu
    //   -logx - R::lbeta(shape3, shape1);
    double logP = log_shape2 + shape3 * logu + shape1 * log1mu
      -logx - lbeta_shape3_shape1;
    return give_log ? logP : std::exp(logP);
  }
  
  
};




// =============================================================================
// Just for reference.
// =============================================================================
double ptrbeta(double q, double shape1, double shape2, double shape3,
               double scale, int lower_tail, int log_p)
{
  if (q <= 0) return 0;
  double logvm, u;
  logvm = shape2 * (std::log(scale) - std::log(q));
  u = std::exp(-log1pexp(logvm));
  if (u > 0.5)
  {
    double u1m = std::exp(-log1pexp(-logvm));
    return R::pbeta(u1m, shape1, shape3, 1 - lower_tail, log_p);
  }
  return R::pbeta(u, shape3, shape1, lower_tail, log_p);
}




struct TrbCdf
{
  double shape1, shape2, shape3, scale;
  double log_scale;
  TrbCdf(){}
  TrbCdf(double a, double b, double c, double d){ reset(a, b, c, d); }
  void reset(double a, double b, double c, double d)
  {
    shape1 = a / c; shape2 = c; shape3 = b / c; scale = d;
    log_scale = std::log(scale);
  }
  
  
  double operator()(double q)
  {
    if (q <= 0) return 0;
    double logvm, u;
    logvm = shape2 * (log_scale - std::log(q));
    u = std::exp(-log1pexp(logvm));
    if (u > 0.5)
    {
      double u1m = std::exp(-log1pexp(-logvm));
      return R::pbeta(u1m, shape1, shape3, 1 - 1, 0);
    }
    return R::pbeta(u, shape3, shape1, 1, 0);
  }
  
};
    












































