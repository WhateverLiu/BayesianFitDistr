// [[Rcpp::export]]
List extractMain(List distlist, NumericVector MDR, 
                 bool normalizeMainPart = false,
                 bool removeZeroPs = false)
{
  List rst(distlist.size());
  NumericVector lm1(distlist.size());
  if (!removeZeroPs)
  {
    for (int i = 0, iend = distlist.size(); i < iend; ++i)
    {
      List d = distlist[i];
      NumericVector val = d[0], P = d[1];
      NumericVector valNew(val.begin() + 1, val.end());
      NumericVector Pnew(P.begin() + 1, P.end());
      // double r = 1.0 / (1.0 - P[0]);
      double r = normalizeMainPart ? 1.0 / (1.0 - P[0]) : 1.0;
      for (auto it = Pnew.begin(); it < Pnew.end(); ++it) *it *= r;
      lm1[i] = MDR[i] * r;
      rst[i] = List::create(Named("val") = valNew, Named("P") = Pnew);
    }
  }
  else
  {
    for (int i = 0, iend = distlist.size(); i < iend; ++i)
    {
      List d = distlist[i];
      NumericVector val = d[0], P = d[1];
      int Nnonzero = 0;
      for (auto it = P.begin() + 1; it < P.end(); ++it) Nnonzero += *it != 0;
      NumericVector valNew(Nnonzero);
      NumericVector Pnew(Nnonzero);
      // double r = 1.0 / (1.0 - P[0]);
      double r = normalizeMainPart ? 1.0 / (1.0 - P[0]) : 1.0;
      int h = 0;
      for (int k = 1, kend = P.size(); k < kend; ++k)
      {
        if (P[k] == 0) continue;
        valNew[h] = val[k];
        Pnew[h] = P[k] * r;
        ++h;
      }
      lm1[i] = MDR[i] * r;
      rst[i] = List::create(Named("val") = valNew, Named("P") = Pnew);
    } 
  }
  
  
  return List::create(Named("conditionalDistrs") = rst, Named("lm1") = lm1);
}



















