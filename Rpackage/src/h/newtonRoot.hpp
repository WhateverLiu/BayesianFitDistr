#pragma once

// Lower bounded by 1e-10.
template <typename Fun> // f.first = value, f.second = gradient.
inline double newtonRoot(int &iter, double x, Fun &f, double eps, int maxit)
{
  double lb = 1e-10, rb = 1e300; // left and right barrier on the varying parameter.
  iter = 0;
  std::pair<double, double> fx(0, 0);
  for (; iter < maxit; ++iter)
  {
    fx = f(x);
    if (fx.first > eps) rb = x;
    else if (fx.first < -eps) lb = x;
    else break;
    double xnew = x - fx.first / fx.second;
    // Rcout << "fx.first = " << fx.first << ", ";
    // Rcout << "fx.second = " << fx.second << ", ";
    // Rcout << "x = " << x << ", xnew = " << xnew << ", lb = " << lb << ", ub = " << rb << "\n";
    if (xnew <= lb) x = (x + lb) * 0.5;
    else if (xnew >= rb) x = (x + rb) * 0.5;
    else x = xnew;
  }
  return x;
}
  

template <typename Fun> // f.first = value, f.second = gradient.
inline double newtonRoot(double x, Fun &f, double eps, int maxit = 100)
{
  int iter;
  return newtonRoot(iter, x, f, eps, maxit);
}
    























