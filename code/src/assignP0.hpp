// [[Rcpp::export]]
void assignP0(NumericVector p, double p0) 
{
  double main = std::accumulate(p.begin() + 1, p.end(), 0.0);
  double r = (1.0 - p0) / main;
  p[0] = p0;
  for (auto it = p.begin() + 1; it < p.end(); ++it) *it *= r;
}


// // [[Rcpp::export]]
// void inplaceNormalize(NumericVector p)
// {
//   double r = 1.0 / std::accumulate(p.begin(), p.end(), 0.0);
//   for (auto it = p.begin(); it != p.end(); ++it) *it *= r;
// }
