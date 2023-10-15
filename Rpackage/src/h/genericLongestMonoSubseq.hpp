#pragma once
#define vec std::vector


// Return the zero-based indices of the longest monotonic subsequence.
template<typename ing, typename num>
struct LongestMonoSubseq
{
  vec<ing> C;
  vec<num> subseq;
  
  
  // Write the indices of the subsequence on v[]. v has to be a space of size
  // at least xend - x. The function returns the pointer that points to the
  // end of the subsequence.
  template <typename Compare>
  ing * run(num*x, num *xend, ing *rst, Compare &&cmp)
  {
    auto xsize = xend - x;
    C.resize(xsize);
    C[0] = 1;
    // Find the size of the longest increasing subsequence ending at x[i].
    // The size is stored in C[i].
    subseq.resize(xsize);
    subseq.resize(1);
    subseq[0] = x[0];
    for (auto i = x-x+1; i < xsize; ++i)
    {
      auto it = std::lower_bound(subseq.begin(), subseq.end(), x[i], cmp);
      if (it >= subseq.end())
      { 
        subseq.emplace_back(x[i]);
        C[i] = subseq.size();
      }
      else
      {
        *it = x[i];
        C[i] = it - subseq.begin() + 1;
      }
    }
    auto whichMax = std::max_element(C.begin(), C.end()) - C.begin();
    rst[0] = whichMax;
    int64_t k = 1;
    for (int64_t i = whichMax - 1; i >= 0; --i)
    { 
      if (C[i] == C[whichMax] - 1)
      {
        rst[k] = i; k += 1;
        whichMax = i;
      }  
    } 
    std::reverse(rst, rst + k);
    return rst + k;
  } 
  
  
  // monoType == 0: increasing subsequence.
  // monoType == 1: nondecreasing subsequence.
  // monoType == 2: decreasing subsequence.
  // monoType == 3: nonincreasing subsequence.
  // Write the indices of the subsequence on rst[]. rst[] MUST be a space of size
  // at least xend - x. The function returns the pointer that points to the
  // end of the subsequence.
  template<int monoType>
  ing * operator()(num *x, num *xend, ing *rst)
  { 
    if (monoType == 0) return run(
      x, xend, rst, [](const num &a, const num &b)->bool { return a < b; });
    if (monoType == 1) return run(
      x, xend, rst, [](const num &a, const num &b)->bool { return a <= b; });
    if (monoType == 2) return run(
      x, xend, rst, [](const num &a, const num &b)->bool { return a > b; });
    if (monoType == 3) return run(
      x, xend, rst, [](const num &a, const num &b)->bool { return a >= b; });
    return nullptr;
  }
};



#undef vec
