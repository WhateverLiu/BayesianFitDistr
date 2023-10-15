#include <random>


// rst[] and [begin, end) shall not overlap.
// rst[] is a buffer of size at least `sampleSize`.
// This could result in 
struct StratifiedSampling
{
  template<typename ing, typename RNG>
  void operator()(ing popuSize, ing sampleSize, ing *rst, 
                bool replace, RNG &rng)
  {
    replace |= sampleSize > popuSize;
    if (!replace and sampleSize == popuSize) 
    {
      std::iota(rst, rst + popuSize, 0);
      return;
    }
    

    double dt = popuSize / double(sampleSize);
    std::uniform_real_distribution<double> U(0, 1);
    
    
    double u = U(rng);
    rst[0] = (0 + u) * dt;
    for (ing i = 1; i < sampleSize; ++i)
    { 
      u = U(rng);
      ing k = (i + u) * dt;
      if (!replace and k == rst[i - 1])
      {
        // Generate a uniform between [ k + 1, (i + 1) * dt ):
        k = ((i + 1) * dt - (k + 1)) * u + (k + 1);
      }
      rst[i] = k;
    } 
  }
};







