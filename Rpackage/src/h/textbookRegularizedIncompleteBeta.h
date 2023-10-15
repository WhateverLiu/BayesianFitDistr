#pragma once

// struct Beta : Gauleg18 {
// Object for incomplete beta function. Gauleg18 provides coefficients for Gauss-Legendre quadrature.
// static const int SWITCH = 3000; When to switch to quadrature method.
// static const double EPS, FPMIN; See end of struct for initializations.

// #define gammln std::lgamma


// constexpr const double BetaEPS = std::numeric_limits<double>::epsilon();
// constexpr const double BetaFPMIN = std::numeric_limits<double>::min() / BetaEPS;


struct NumericalRecipesRIB 
{
  
  
  static double betacf(const double a, const double b, const double x, 
                       const double fpmin, const double EPS) 
  {
    int m, m2;
    double aa, c, d, del, h, qab, qam, qap;
    qab = a + b; // These q’s will be used in factors that
    qap = a + 1.0; // occur in the coefficients (6.4.6).
    qam = a - 1.0;
    c = 1.0; // First step of Lentz’s method.
    d = 1.0 - qab * x / qap;
    if (std::abs(d) < fpmin) d = fpmin;
    d = 1.0 / d;
    h = d;
    for (m = 1; m < 10000; ++m) 
    {
      m2 = 2 * m;
      aa = m * (b - m) * x / ( (qam + m2) * (a + m2) );
      d = 1.0 + aa * d; // One step (the even one) of the recurrence.
      if (std::abs(d) < fpmin) d = fpmin;
      c = 1.0 + aa / c;
      if (std::abs(c) < fpmin) c = fpmin;
      d = 1.0 / d;
      h *= d * c;
      aa = -(a + m) * (qab + m) * x / ( (a + m2) * (qap + m2) );
      d = 1.0 + aa * d; // Next step of the recurrence (the odd one).
      if (std::abs(d) < fpmin) d = fpmin;
      c = 1.0 + aa / c;
      if (std::abs(c) < fpmin) c = fpmin;
      d = 1.0 / d;
      del = d * c;
      h *= del;
      if (std::abs(del - 1.0) <= EPS) break; // Are we done?
    }
    return h;
  }
  
  
  /*
   double betaiapprox(double a, double b, double x) 
   {
   // Incomplete beta by quadrature. Returns Ix.a; b/. User should not call directly.
   int j;
   double xu, t, sum, ans;
   double a1 = a - 1.0, b1 = b - 1.0, mu = a / (a + b);
   double lnmu = std::log(mu), lnmuc = std::log(1.0 - mu);
   t = std::sqrt(a * b / ( (a + b) * (a + b) * (a + b + 1.0)) );
   if (x > a / (a + b)) // Set how far to integrate into the tail:
   { 
   if (x >= 1.0) return 1.0;
   xu = std::min(1.0, std::max( mu + 10.0 * t, x + 5.0 * t ));
   } 
   else 
   {
   if (x <= 0.0) return 0.0;
   xu = std::max(0.0, std::min(mu - 10.0 * t, x - 5.0 * t));
   }
   sum = 0;
   
   
   double y[18] = {
   0.0021695375159141994,
   0.011413521097787704,0.027972308950302116,0.051727015600492421,
   0.082502225484340941, 0.12007019910960293,0.16415283300752470,
   0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
   0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
   0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
   0.87126389619061517, 0.95698180152629142};
   
   
   double w[18] = {
   0.0055657196642445571,
   0.012915947284065419,0.020181515297735382,0.027298621498568734,
   0.034213810770299537,0.040875750923643261,0.047235083490265582,
   0.053244713977759692,0.058860144245324798,0.064039797355015485,
   0.068745323835736408,0.072941885005653087,0.076598410645870640,
   0.079687828912071670,0.082187266704339706,0.084078218979661945,
   0.085346685739338721,0.085983275670394821};
   
   
   for (j = 0; j < 18; j++) // Gauss-Legendre.
   { 
   t = x + (xu - x) * y[j];
   sum += w[j] * std::exp( a1 * (std::log(t) - lnmu) + 
   b1 * (std::log(1 - t) - lnmuc) );
   }
   ans = sum * (xu - x) * std::exp( a1 * lnmu - std::lgamma(a) + 
   b1 * lnmuc - std::lgamma(b) + std::lgamma(a + b));
   return ans > 0.0 ? 1.0 - ans : -ans;
   }
   */
  
  
  static double betai(const double x, const double a, const double b,
                      const double fpmin, const double eps) 
  {
    double bt;
    if (x == 0.0 || x == 1.0) return x;
    // if (a > 3000 && b > 3000) return betaiapprox(a, b, x);
    bt = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) + 
      a * std::log(x) + b * std::log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0) ) 
      return bt * betacf(a, b, x, fpmin, eps) / a;
    return 1.0 - bt * betacf(b, a, 1.0 - x, fpmin, eps) / b;
  }
  
  
  
  
  static double invbetai(double p, double a, double b, 
                         const double fpmin, const double eps) 
  {
    // Inverse of incomplete beta function. Returns x such that Ix.a; b/ D p for argument p
    // between 0 and 1.
    double pp, t, u, err, x, al, h, w, afac, a1 = a - 1.0, b1 = b - 1.0;
    int j;
    if (p <= 0.0) return 0.0;
    else if (p >= 1.0) return 1.0;
    else if (a >= 1.0 && b >= 1.0) // Set initial guess. See text.
    { 
      pp = p < 0.5 ? p : 1.0 - p;
      t = std::sqrt(-2.0 * std::log(pp));
      x = (2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t;
      if (p < 0.5) x = -x;
      // al = (SQR(x) - 3.0) / 6.0;
      al = (x * x - 3.0) / 6.0;
      h = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0) );
      w = (x * std::sqrt( al + h) / h) - (1.0 / (2.0 * b - 1) - 1.0 / 
        (2.0 * a - 1.0)) * (al + 5.0 / 6.0 - 2.0 / (3.0 * h));
      x = a / ( a + b * std::exp(2.0 * w));
    } 
    else 
    {
      double lna = std::log( a / (a + b) ), lnb = std::log( b / (a + b) );
      t = std::exp(a * lna) / a;
      u = std::exp(b * lnb) / b;
      w = t + u;
      if (p <  t / w) x = std::pow( a * w * p, 1.0 / a);
      else x = 1.0 - std::pow(b * w * (1.0 - p), 1.0 / b);
    }
    
    
    afac = -std::lgamma(a) - std::lgamma(b) + std::lgamma(a + b);
    
    
    Rcout << "in invbeta, x = " << x << "\n";
    
    
    // const double EPS = 1.0e-8;
    for (j = 0; j < 100; ++j) 
    {
      if (x == 0.0 || x == 1.0) return x; // a or b too small for accurate calculation.
      err = betai(a, b, x, fpmin, eps) - p;
      t = std::exp( a1 * std::log(x) + b1 * std::log(1.0 - x) + afac ); 
      u = err / t; // Halley:
      x -= (t = u / (1.0 - 0.5 * std::min(1.0, u * (a1 / x - b1 / (1.0 - x)))));
      if (x <= 0.0) x = 0.5 * (x + t); // Bisect if x tries to go neg or > 1.
      if (x >= 1.0) x = 0.5 * (x + t + 1.0);
      if (std::abs(t) < 1e-8 * x and j > 0) break;
    }
    Rcout << "in invbeta, returned x = " << x << "\n";
    Rcout << "in invbeta, number of iterations = " << j << "\n";
    return x;
  } 
  // };
  
  
};



// #undef gammln






















