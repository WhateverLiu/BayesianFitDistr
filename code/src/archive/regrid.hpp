// #include <Rcpp.h>
// using namespace Rcpp;


#define losstype double
#define probtype double


void normalize(double *p, double *pend)
{
  double *presv = p;
  double sum = 0;
  for (; p < pend; ++p)
  {
    *p = std::max(0.0, *p);
    sum += *p;
  }
  double r = 1.0 / sum;
  p = presv;
  for (; p < pend; ++p) *p *= r;
}


bool rg4_(
    losstype *inLoss, probtype *inProb, int inSize, losstype *outLoss,
    probtype *outProb, int outSize, double *safe) 
{
  
  bool unableToWork = 0;
  std::copy(outProb, outProb + outSize, safe);
  
  
  if (outSize < 5 or inSize < 5) 
  {
    unableToWork = 1;
  } 
  else
  {
    
    int n = 0;
    losstype x0 = outLoss[0];
    losstype x3 = outLoss[outSize - 1];
    probtype p0 = 0;
    probtype p1 = 0;
    probtype p2 = 0;
    probtype p3 = 0;
    probtype p = 0;
    losstype x = 0;
    losstype x1 = 0;
    losstype x2 = 0;
    for (int i = 0, iend = outSize - 1; i < iend; ) 
    {
      if (std::abs(inLoss[n] - outLoss[i]) < 1e-10) {
        outProb[i] += inProb[n];
        Rcout << "???" << inProb[n] << "\n";
        ++n;
        if (n >= inSize) {
          break;
        }
      } else {
        x1 = outLoss[i];
        x2 = outLoss[i + 1];
        while (n < inSize and (x2 - inLoss[n]) > 1e-10)
        {
          x = inLoss[n];
          p = inProb[n];
          double tmp = p * (x3 - x) * (x0 - x) /
            ((x3 - x1) * (x0 - x1) * (x2 - x) +
              (x3 - x2) * (x0 - x2) * (x - x1));
          p1 = tmp * (x2 - x);
          p2 = tmp * (x - x1);
          p0 = (p * (x3 - x) - (x3 - x1) * p1 - (x3 - x2) * p2) / (x3 - x0);
          p3 = p - p0 - p1 - p2;
          ++n;
          outProb[0] += p0;
          outProb[i] += p1;
          outProb[i + 1] += p2;
          outProb[outSize - 1] += p3;
        }
        i++;
      }
    }
    
    
    if (std::abs(outLoss[outSize - 1] - inLoss[inSize - 1]) < 1e-10)
      outProb[outSize - 1] += inProb[inSize - 1];
    
    
    if (outProb[0] >= 0 and outProb[outSize - 1] >= 0)
    {
      normalize(outProb, outProb + outSize);
      return 0; 
    }
    
    
    int i = 0;
    int j = outSize - 1;
    
    
    while (true) 
    {
      if (outProb[i] < 0) {
        x = outLoss[i];
        x1 = outLoss[i + 1];
        x2 = outLoss[i + 2];
        x3 = outLoss[j];
        
        p = outProb[i];
        
        p1 = (x - x3) * (x - x2) / ((x1 - x3) * (x1 - x2)) * p;
        p2 = (x - x1) * (x - x3) / ((x2 - x3) * (x2 - x1)) * p;
        p3 = p - p1 - p2;
        outProb[i + 1] += p1;
        outProb[i + 2] += p2;
        outProb[j] += p3;
        outProb[i] = 0;
        ++i;
      }
      if (outProb[i] >= 0 and outProb[j] >= 0) {
        break;
      }
      
      if (j - i < 3) {
        unableToWork = true;
        break;
      }
      
      if (outProb[j] < 0) {
        x = outLoss[j];
        x1 = outLoss[j - 1];
        x2 = outLoss[j - 2];
        x3 = outLoss[i];
        p = outProb[j];
        
        p1 = (x - x3) * (x - x2) / ((x1 - x3) * (x1 - x2)) * p;
        p2 = (x - x1) * (x - x3) / ((x2 - x3) * (x2 - x1)) * p;
        p3 = p - p1 - p2;
        outProb[j - 1] += p1;
        outProb[j - 2] += p2;
        outProb[i] += p3;
        outProb[j] = 0;
        --j;
      }
      if (outProb[i] >= 0 and outProb[j] >= 0) {
        break;
      }
      if (j - i < 3) {
        unableToWork = true;
        break;
      }
    }
  }
  
  
  if (unableToWork) 
  {
    int j = 0;
    std::copy(safe, safe + outSize, outProb);
    double tmp, interval;
    losstype *ng_st = outLoss, *x_st = inLoss;
    probtype *ngp_st = outProb, *xp_st = inProb;
    int ngsize = outSize, sizex = inSize;
    for (int i = 1; i < ngsize; ++i) 
    {
      interval = ng_st[i] - ng_st[i - 1];
      while (j < sizex and x_st[j] <= ng_st[i]) 
      {
        tmp = xp_st[j] / interval * ( ng_st[i] - x_st[j] );
        ngp_st[i - 1] += tmp;
        ngp_st[i] += xp_st[j] - tmp;
        ++j;
      }
    }
  }
  
  
  normalize(outProb, outProb + outSize);
  return unableToWork;
}


#undef losstype
#undef probtype


// Will match the results from keyALGs.
// [[Rcpp::export]]
DataFrame rg4(DataFrame X, NumericVector ngd) 
{
  NumericVector val = X[0], P = X[1], newgridP(ngd.size());
  vec<double> safe(newgridP.size());
  rg4_(&val[0], &P[0], val.size(), &ngd[0], &newgridP[0], ngd.size(), &safe[0]); 
  return DataFrame::create(Named("val") = ngd, Named("P") = newgridP);
}


// [[Rcpp::export]]
DataFrame empc(NumericVector x, int n)
{
  n = std::max(n, 2);
  int size = x.size();
  if (size <= 1) return DataFrame::create(Named("val") = 0.0, Named("P") = 1.0);
  vec<double> cnt;
  bool sorted = true;
  for (int i = 1; i < size; ++i)
  {
    if (x[i] < x[i - 1]) { sorted = false; break;  }
  }
  double *y = nullptr, *p = nullptr;
  if (sorted) 
  {
    cnt.resize(size * 2);
    y = &x[0];
    p = &cnt[0];
  }
  else 
  {
    cnt.resize(size * 3);
    y = &cnt[0];
    std::copy(x.begin(), x.end(), y);
    std::sort(y, y + size);
    p = y + size;
  }
  double *safe = p + size;
  double ymin = y[0], ymax = y[size - 1];
  double delta = (ymax - ymin) / (n - 1);
  NumericVector ngd(n);
  for (int i = 0; i < size; ++i) p[i] = 1.0 / size;
  for (int i = 0, iend = n - 1; i < iend; ++i) ngd[i] = ymin + delta * i;
  ngd[n - 1] = ymax;
  NumericVector pnew(n);
  rg4_(y, p, size, &ngd[0], &pnew[0], ngd.size(), safe);
  return DataFrame::create(Named("val") = ngd, Named("P") = pnew);
}


// [[Rcpp::export]]
DataFrame makeHist(NumericVector x, int n)
{
  double ymin = 1e300, ymax = -1e300;
  for (int i = 0, iend = x.size(); i < iend; ++i)
  {
    ymin = std::min(ymin, x[i]);
    ymax = std::max(ymax, x[i]);
  }
  double dy = (ymax - ymin) / (n - 1);
  NumericVector sp(n), p(n);
  for (int i = 0; i < n; ++i) sp[i] = ymin + i * dy;
  ymin -= dy / 2;
  double r = 1.0 / x.size();
  for (int i = 0, iend = x.size(); i < iend; ++i)
  {
    int k = (x[i] - ymin) / dy;
    p[k] += r;
  }
  return DataFrame::create(Named("val") = sp, Named("P") = p);
}



























// explain option II
// change notations
// counter examples




