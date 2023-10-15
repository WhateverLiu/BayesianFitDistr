#include <Rcpp.h>
using namespace Rcpp;


template<typename num>
struct A
{
  template<int y>
  num run(num x) { return x + y; }
};


template<typename num>
struct B
{
  A<num> a;
  // A<double> a;
  num run(num x)
  // double run(double x)
  {
    if (x > 1) return a.template run<1> (x);
    return a.template run<0> (x);
  }
};


// // [[Rcpp::export]]
double test(double x)
{
  B<double> b;
  // B b;
  return b.run(x);
}




