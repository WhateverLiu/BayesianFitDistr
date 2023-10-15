#define vec std::vector


// C++ standard requires allocation functions (malloc() and operator new()) to 
// allocate memory suitably aligned for any standard type. As these functions 
// don't receive the alignment requirement as an argument, in practice it means 
// that the alignment for all allocations is the same, and is that of a standard 
// type with the largest alignment requirement.


// Every thread should have its own VecPool.
class VecPool
{
  vec<vec<char> > X; // Initialization is std::size_t byte aligned.
  
  
  template<typename T>
  vec<char> T2char(vec<T> &x)
  { 
    static_assert(sizeof(x) == sizeof(vec<char>));
    x.clear(); // Call destructors on all elements but maintain capacity.
    vec<char> rst;
    std::size_t mem[3]; // Start swapping.
    mmcp(mem,  &rst, sizeof(rst));
    mmcp(&rst, &x,   sizeof(x));
    mmcp(&x,   mem,  sizeof(x));
    return rst;
  }
  
  
  template<typename T>
  vec<T> char2T(vec<char> &x)
  { 
    static_assert(sizeof(x) == sizeof(vec<T>));
    x.clear();
    vec<T> rst;
    std::size_t mem[3];
    mmcp(mem, &x, sizeof(x));
    // Adjust the capacity pointer to ensure the capacity is a multiple 
    //   of sizeof(T).
    mem[2] -= (mem[2] - mem[0]) % sizeof(T);
    mmcp(&x, &rst, sizeof(rst));
    mmcp(&rst, mem, sizeof(rst));
    return rst;
  }
  
  
public:
  
  
  template<typename DestPtr, typename SourPtr>
  void mmcp(DestPtr dest, SourPtr src, std::size_t Nbyte) 
  { 
    std::memcpy((void*)dest, (void*)src, Nbyte); 
  }
  
  
  vec<vec<char> > &getX() { return X; }
  
  
  VecPool()
  {
    static_assert( sizeof(std::size_t) == sizeof(void*) );
    
    
    std::srand(std::time(nullptr));
    
    
    // Test if std::vector is implemented in the way that can be exploited.
    if (true)
    {
      int fullSize = std::rand() % 31 + 7;
      int subSize = std::rand() % (fullSize / 2) + 1;
      
      
      typedef std::tuple<char, char, char> tupe; // Arbitrarily selected type.
      vec<tupe> a(fullSize); // Just for creating a vector of an arbitrary type.
      a.resize(subSize); // Downsize the vector.
      
      
      // Vector header must contain exactly 3 pointers, e.g. p1, p2, p3.
      static_assert( sizeof(a) == sizeof(std::size_t) * 3 );
      std::size_t ptr[3];
      mmcp(ptr, &a, sizeof(a));
      
      
      // The 3 pointers also must satisfy the following constraints:
      if ( !( 
          ptr[0] + sizeof(tupe) * subSize  == ptr[1] and
             ptr[0] + sizeof(tupe) * fullSize == ptr[2]
      )) 
        throw std::runtime_error(
            "Layout of std::vector header is not supported.");
    }
    
    
    
    // We check a few things. First, after we returned a container to VecPool,
    // does VecPool actuall have the same container. 
    // VecPool and returned it to
    // float arr[3];
    // for (int i = 0; i < 3; ++i) arr[i] = std::rand() * 3.14159;
    auto a = give<float>();
    a.resize(5);
    a.resize(3);
    std::copy(arr, arr + 3, a.begin());
    void *aptr = &a[0];
    
    
    // std::tuple<short, short, short> brr[5];
    // mmcp(brr, arr, sizeof(float) * 3);
    auto b = give<std::tuple<short, short, short> >();
    b.resize(5);
    // std::copy(brr, brr + 5, b.begin());
    void *bptr = &b[0];
    
    
    this->recall(b);
    this->recall(a);
    
    
    // if (std::memcmp(arr, &X[1][0], sizeof(float) * 3) != 0)
    //   throw std::runtime_error("VecPool is not correctly built --- 0.");
    // // else Rcout << "0 yes!\n"; // DO NOT DELETE COMMENT!
    // 
    // 
    // if (std::memcmp(brr, &X[0][0], sizeof(
    //     std::tuple<short, short, short>) * 5) != 0)
    //   throw std::runtime_error("VecPool is not correctly built --- 1.");
    // // else Rcout << "1 yes!\n"; // DO NOT DELETE COMMENT!
    
    
    if ( (void*)(&X[1][0]) != aptr )
      throw std::runtime_error("VecPool is not correctly built --- 2.");
    // else Rcout << "2 yes!\n"; // DO NOT DELETE COMMENT!
    
    
    if ( (void*)(&X[0][0]) != bptr )
      throw std::runtime_error("VecPool is not correctly built --- 3.");
    // else Rcout << "3 yes!\n"; // DO NOT DELETE COMMENT!
    
    
    auto c = give<std::tuple<char, char, char> >();
    c.resize(2);
    constexpr const int byteSize = sizeof(std::tuple<char, char, char>) * 2;
    void *cptr = &c[0];
    double num = 1.0 / 3.14159265;
    mmcp(&c[0], &num, byteSize);
    char crr[byteSize];
    mmcp(crr, &num, byteSize);
    recall(c);
    if ( (void*)(&X.back().front()) != cptr )
      throw std::runtime_error("VecPool is not correctly built --- 4.");
    // else Rcout << "4 yes!\n"; // DO NOT DELETE COMMENT!
    
    
    if ( std::memcmp(crr, &X.back().front(), byteSize) != 0)
      throw std::runtime_error("VecPool is not correctly built --- 5.");
    // else Rcout << "5 yes!\n"; // DO NOT DELETE COMMENT!
  }
  
  
  template<typename T>
  vec<T> give()
  {
    if (X.size() == 0) return vec<T>(0);
    vec<char> rst;
    rst.swap(X.back());
    rst.resize(0);
    X.pop_back();
    return char2T<T> (rst);
  }
  
  
  template<typename T>
  void recall(vec<T> &x) // invalidates everything related to x.
  {
    X.emplace_back(T2char(x));
  } 
  
  
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
  
};


#undef vec




// [[Rcpp::export]]
void testVecPool() {
  VecPool vp;
  Rcout << vp.getX().back().size() << ", " << vp.getX().back().capacity() << "\n";
  auto x = vp.give<double> ();
  Rcout << x.size() << ", " << x.capacity() << "\n";
  x.emplace_back(1.1);
  x.emplace_back(1.2);
  x.emplace_back(1.3);
  for (int i = 0, iend = x.size(); i < iend; ++i)
    Rcout << x[i] << ", ";
  Rcout << "\n";
  std::size_t p[3];
  std::memcpy((void*)p, (void*)(&x), sizeof(x));
  Rcout << p[0] << ", " << p[1] << ", " << p[2] << "\n";
  vp.recall(x);
  Rcout << vp.getX().back().size() << ", " << vp.getX().back().capacity() << "\n";
}






