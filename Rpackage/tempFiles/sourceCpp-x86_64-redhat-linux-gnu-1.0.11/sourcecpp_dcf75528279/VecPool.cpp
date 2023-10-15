// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;


#define vec std::vector


template<typename DestPtr, typename SourPtr>
inline void mmcp(DestPtr dest, SourPtr src, std::size_t Nbyte) 
{ 
  std::memcpy((void*)dest, (void*)src, Nbyte); 
}


// =============================================================================
// Swap two std::vectors of different types. After swapping, both vectors'
// sizes will become 0. For each of them, the capacity pointer could be shifted
// to ensure that the capacity is a multiple of sizeof(type). 
// No allocation will occur. The vectors maintain their original element types.
// =============================================================================
struct DarkSwapVec
{
  DarkSwapVec()
  {
    static_assert( sizeof(std::size_t) == sizeof(void*) );
    
    
    // Test if std::vector is implemented in the way that can be exploited.
    constexpr const unsigned fullSize = 31;
    constexpr const unsigned subSize  = 13;
    
    
    typedef std::tuple<char, char, char> tupe; // Arbitrarily selected type.
    vec<tupe> a(fullSize); // Just for creating a vector of an arbitrary type.
    a.resize(subSize); // Downsize the vector.
    
    
    // Vector header must contain exactly 3 pointers, e.g. p1, p2, p3.
    static_assert( sizeof(a) == sizeof(std::size_t) * 3 );
    auto ptr = (std::size_t*)(&a);
    bool sizePointerCool     = ptr[0] + sizeof(tupe) * subSize  == ptr[1];
    bool capacityPointerCool = ptr[0] + sizeof(tupe) * fullSize == ptr[2];
    
    
    if (!sizePointerCool or !capacityPointerCool)
      throw std::runtime_error("Layout of std::vector header is unsupported.");
  }
  
  
  template<typename S, typename T>
  void operator()(vec<S> &x, vec<T> &y)
  {
    static_assert( sizeof(std::size_t) == sizeof(void*) );
    static_assert( sizeof(x) == sizeof(y) );
    static_assert( sizeof(x) == sizeof(std::size_t) * 3 );
    x.clear(); // Enforce destructor on elements if type is nontrivial.
    y.clear();
    // Rcout << "x.capacity() = " << x.capacity() << "\n";
    // Rcout << "y.capacity() = " << y.capacity() << "\n";
    
    
    std::size_t xt[3], yt[3];
    mmcp(xt, &x, sizeof(x));
    mmcp(yt, &y, sizeof(y));
    
    
    // for (int i = 0;i < 3; ++i) Rcout << xt[i] << " ";
    // Rcout << "\n";
    // for (int i = 0;i < 3; ++i) Rcout << yt[i] << " ";
    // Rcout << "\n\n";
    
    
    xt[2] -= (xt[2] - xt[0]) % sizeof(T);
    yt[2] -= (yt[2] - yt[0]) % sizeof(S);
    
    
    // for (int i = 0;i < 3; ++i) Rcout << xt[i] << " ";
    // Rcout << "\n";
    // for (int i = 0;i < 3; ++i) Rcout << yt[i] << " ";
    // Rcout << "\n\n";
    
    
    mmcp(&x, yt, sizeof(x));
    mmcp(&y, xt, sizeof(y));
    
    
    // Rcout << "x.capacity() = " << x.capacity() << "\n";
    // Rcout << "y.capacity() = " << y.capacity() << "\n";
    // Rcout << "\n\n";
  }
};




template <bool check = false>
struct VecPool
{
  DarkSwapVec darkSwap;
  vec<vec<char> > X;
  
  
  void examineExample()
  {
    auto a = give<float>();
    a.resize(3);
    void *aptr = a.data();
    
    
    auto b = give<std::tuple<short, short, short> >();
    b.resize(5);
    void *bptr = b.data();
    
    
    this->recall(b);
    this->recall(a);
    
    
    if ( (void*)(X.back().data()) != aptr )
      throw std::runtime_error("Last in pool is a different container.");
    
    
    if ( (void*)(X[X.size() - 2].data()) != bptr )
      throw std::runtime_error("Second last in pool is a different container.");
    
    
    typedef std::tuple<char, char, char, char, char, char, char> tupe;
    auto c = give<tupe>();
    if (c.capacity() % sizeof(tupe) != 0) throw std::runtime_error(
      "Container capacity is not a multiple of element byte size.");
    c.resize(2);
    auto *cptr = c.data();
    recall(c);
    if ( (void*)(X.back().data()) != cptr )
      throw std::runtime_error("Last in pool should not change but did.");
  }
  
  
  // ===========================================================================
  // Dispatch and recall containers.
  // ===========================================================================
  template<typename T>
  vec<T> give()
  {
    if (X.size() == 0) return vec<T> (0);
    vec<T> rst;
    if constexpr (check) 
      std::cout << "VecPool::give() X.back().capacity() = " + 
        std::to_string(X.back().capacity()) + "\n";
    darkSwap(X.back(), rst);
    if constexpr (check)
      std::cout << "VecPool::give() output.capacity() = " + 
        std::to_string(rst.capacity()) + ", type size = " + 
        std::to_string(sizeof(T)) + "\n\n";
    X.pop_back();
    return rst;
  } 
  
  
  template<typename T>
  void recall(vec<T> &x) // invalidates everything related to x.
  { 
    x.clear();
    if constexpr (check) 
      std::cout << "VecPool::recall() input.capacity() = " + 
        std::to_string(x.capacity()) + ", type size = " + 
        std::to_string(sizeof(T)) + "\n";
    vec<char> rst;
    darkSwap(rst, x);
    // =========================================================================
    // Do not directly do X.emplace_back(rst). This will have no effect since
    //   rst.size() == 0. emplace_back ignores zero-size vectors. Tricky!
    // =========================================================================
    X.emplace_back(vec<char>(0));
    X.back().swap(rst);
    if constexpr (check) 
      std::cout << "VecPool::recall() X.back().capacity() = " + 
        std::to_string(X.back().capacity()) + "\n\n";
  }
};




// ===========================================================================
// Code pattern goes like this:
// ===========================================================================
// VecPool vp;
// auto a = vp.give<int> ();
// auto b = vp.give<double>();
// auto c = vp.give<std::tuple<short, short, short> >();
// 
// ......
// 
// vp.recall(c); // Notice the reversed order of recall. Not necessary but good to do.
// vp.recall(b);
// vp.recall(a);
// ===========================================================================


// [[Rcpp::export]]
int test(int k)
{
  VecPool<true> vp;
  vec<vec<double> > a(20);
  for (int i = 0, iend = a.size(); i < iend; ++i)
    a[i].resize(i + 7);
  vp.recall(a);
  struct Tmp { char a[7]; 
    Tmp() { for (int i = 0; i < 7; ++i) a[i] = i;  } };
  auto b = vp.give<Tmp>();
  for (int i = 0; i < k; ++i)
    b.emplace_back( Tmp() );
  std::cout << "b.capacity() = " + std::to_string(b.capacity()) + "\n";
  return b.back().a[3];
}




#undef vec






#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// test
int test(int k);
RcppExport SEXP sourceCpp_35_test(SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(test(k));
    return rcpp_result_gen;
END_RCPP
}
