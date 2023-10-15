

// valnew and Pnew have the same sizes as the old PMF's.
template <typename RNG>
void samplePmf(double *val, double *P, int size, 
               double *valnew, double *Pnew, int &sizenew,
               int Nsample, RNG &rng)
{
  double delta = 1.0 / Nsample;
  std::vector<double> cdf(Nsample);
  std::uniform_real_distribution<double> U(0.0, delta);
  for (int i = 0; i < Nsample; ++i) cdf[i] = delta * i + U(rng);
  double psum = 0;
  for (int i = 0; i < size; ++i)
  {
    
  }
    
    
}



// [[Rcpp::export]]
NumericVector asd(NumericVector x) {
  return x * 2;
}






