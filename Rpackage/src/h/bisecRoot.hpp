#pragma once

// Assume f is increasing, and the domin is either nonpositive or nonnegative.
template <typename Fun>
double bisecRoot(double x, Fun &f, double eps, int maxit = 100) // x is initial value.
{
  double fval = f(x);
  if (std::abs(fval) < eps) return x;
  double left, right;
  int it = 0;
  if (fval < 0)
  {
    x *= 2;
    while ( it < maxit and f(x) < 0) { x *= 2; ++it; }
    left = x * 0.5;
    right = x;
  }
  else
  {
    x *= 0.5;
    while ( it < maxit and f(x) > 0 ) { x *= 0.5; ++it; }
    left = x;
    right = x * 2;
  }
  
  
  double mid = 0, fmid;
  for (; it < maxit; ++it)
  {
    mid = (left + right) * 0.5;
    fmid = f(mid);
    if (fmid <= -eps) left = mid;
    else if (fmid >= eps) right = mid;
    else break;    
  }
  // Rcout << fmid << ", ";
  // Rcout << "Iteration = " << it << "\n";
  
  
  return mid;
}











