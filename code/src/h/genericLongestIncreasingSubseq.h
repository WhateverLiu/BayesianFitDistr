

// Write the longest increasing subsequence's indices into rst[], and 
//   return its size.
// rst[] is a buffer of size at least end - begin.
// Stype is the type of the sequence element.


template <typename Stype>
struct LongestIncreasingSubseq
{
  std::vector<Stype> subseq;
  template <typename ing, typename Value, typename LessThanValue, 
            bool strictlyIncreasing = true>
  ing operator()(Stype *begin, Stype *end, const Value &val, 
               LessThan lessthan, ing *rst)
  {
    
  }
};
ing longestIncreasingSubseq(Stype *begin, Stype *end, const Value &val, 
                            LessThan lessthan, ing *rst)
{

  
  
  subseq[0] = nodevec[0];
  subseq.resize(1);
  auto &x = nodevec;
  for (int i = 1, iend = x.size(); i < iend; ++i)
  {
    auto it = std::lower_bound(subseq.begin(), subseq.end(), x[i], lessthan);
    if (it == subseq.end())
    {
      subseq.push_back(x[i]);
      C[i] = subseq.size();
    }
    else
    {
      *it = x[i];
      C[i] = it - subseq.begin() + 1;
    }
  }
  int whichMax = std::max_element(C, C + x.size()) - C;
  int *v = &indicesWanted.back();
  int *vend = v + 1;
  *v = whichMax;
  --v;
  for (int i = whichMax - 1; i >= 0; --i)
  {
    if (C[i] == C[whichMax] - 1)
    {
      *v = i;
      whichMax = i;
      --v;
    }
  }
  longestMonoSeq = v;
  longestMonoSeqEnd = vend;
  
  

}



