// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../h/darkSwapVec.hpp"


// [[Rcpp::export]]
void test() {
  DarkSwapVec dsv;
  std::vector<char> a(11);
  std::vector<double> b(13);
  dsv(a, b);
  
  
  // Rcout << ((std::size_t*)(&a))[0] << ", ";
  // Rcout << ((std::size_t*)(&a))[1] << ", ";
  // Rcout << ((std::size_t*)(&a))[2] << ", ";
  // Rcout << "\n\n";
  // 
  // 
  // Rcout << ((std::size_t*)(&b))[0] << ", ";
  // Rcout << ((std::size_t*)(&b))[1] << ", ";
  // Rcout << ((std::size_t*)(&b))[2] << ", ";
  // Rcout << "\n\n";
  
  
  Rcout << "\n";
  Rcout << a.size() << ", " << a.capacity() << "\n";
  Rcout << b.size() << ", " << b.capacity() << "\n";
}







