#pragma once

struct PlainGD
{
  
  template <typename objFun>
  std::pair<double, int> operator()(
      double *x, // Result parameters will be written into the initialization vector.
      int dim, 
      objFun &f, 
      double learningRate = 1e-4, // learning rate.
      int maxit = 1000,
      double *LB = nullptr, 
      double *UB = nullptr,
      double hgrad = 0, // If zero, will be adjusted inside function.
      bool centralDiff = true,
      double eps = 1e-5, 
      double relaEps = 1e-5
  )
  {
    if (centralDiff) return run<objFun, true>(
        x,
        dim, 
        f, 
        learningRate,
        maxit,
        LB,
        UB,
        hgrad,
        eps,
        relaEps
    );
    
    
    return run<objFun, false>(
        x,
        dim, 
        f, 
        learningRate,
        maxit,
        LB,
        UB,
        hgrad,
        eps,
        relaEps
    );
  }
  
  
  template <typename objFun, bool centralDiff>
  std::pair<double, int> run (
      double *x, // Result parameters will be written into the initialization vector.
      int dim, 
      objFun &f, 
      double learningRate = 1e-4,
      int maxit = 1000,
      double *LB = nullptr, 
      double *UB = nullptr,
      double hgrad = 0, // If zero, will be adjusted inside function.
      double eps = 1e-5, 
      double relaEps = 1e-5
  // =======================================================================
  )
  {
    
    // Amount recommended in paper. MachinePrecision ^ (1/3) and 
    //   MachinePrecision ^ (1/2).
    if (centralDiff) hgrad = hgrad == 0 ? 6.05545445239334e-06 : hgrad;
    else hgrad = hgrad == 0 ? 1.49011611938477e-08 : hgrad;
    
    
    if (LB == nullptr) { lb.assign(dim, -1e300); LB = lb.data(); }
    if (UB == nullptr) { ub.assign(dim,  1e300); UB = ub.data(); }
    grad.resize(dim);
    double *g = grad.data();
    double fval = f(x, g, dim);
    
    
    if (!std::isfinite(fval)) throw std::invalid_argument(
        "Bad initialization. Function value not finite.");
    
    
    int it = 1;
    // double lrResv = learningRate;
    // Rcout << "maxit = " << maxit << "\n";
    for (; it < maxit; ++it)
    {
      
      
      // Rcout << fval << ", " << learningRate << "\n";
      
      
      if (g == nullptr) // Compute gradient via finite differencing.
      {
        for (int i = 0; i < dim; ++i) // Finite difference for gradient.
        { 
          double safe = x[i];
          double dx = std::abs(x[i]) * hgrad;
          x[i] += dx;
          double fhigh = f(x, g, dim);
          if (!centralDiff) grad[i] = (fhigh - fval) / dx;
          else
          {  
            x[i] = safe - dx;
            double flow = f(x, g, dim);
            grad[i] = (fhigh - flow) / (2 * dx);
          }
          x[i] = safe;
        }
      }
      
      
      // Rcout << "gradient = ";
      // for (int k = 0; k < dim; ++k) Rcout << grad[k] << ", ";
      // Rcout << std::endl << std::endl;

      
      double gnorm = std::sqrt(std::inner_product(
        grad.begin(), grad.end(), grad.begin(), 0.0));
      double xnorm = std::sqrt(std::inner_product(x, x + dim, x, 0.0));
      if ( gnorm < xnorm * relaEps or gnorm < eps ) break;
      
      
      xresv.assign(x, x + dim);
      gradTmp.resize(dim);
      // if (g != nullptr) gradResv.assign(grad.begin(), grad.end());
      // gradResv.assign(grad.begin(), grad.end());
      
      
      while (true) // Update attempt.
      // for (int infIter = 0; infIter < 10; ++infIter)
      {
        // Rcout << "learning rate = " << learningRate << ", ";
        // Rcout << "fval = " << fval << std::endl;
        // Rcout << "gradient = ";
        // for (int k = 0; k < dim; ++k) Rcout << grad[k] << ", ";
        
        
        for (int i = 0; i < dim; ++i)
        {
          x[i] = std::max(LB[i], std::min(
            UB[i], x[i] - grad[i] * learningRate));
        }
        g = gradTmp.data();
        double fvalnew = f(x, g, dim);
        // Rcout  << "fvalnew = " << fvalnew << std::endl; 
        if (!std::isfinite(fvalnew) or fvalnew > fval)
        {
          std::copy(xresv.begin(), xresv.end(), x);
          // if (g != nullptr) gradResv.swap(grad);
          // grad.assign(gradResv.begin(), gradResv.end());
          learningRate *= 0.5;
        }
        else 
        {
          // learningRate = lrResv;
          fval = fvalnew;
          learningRate *= 2;
          grad.swap(gradTmp);
          break;
        }
      }
      
      
      // Rcout << "\n\n\n";
      
    }
    
    
    return std::pair<double, int> (fval, it);
  }
  
  
  std::vector<double> lb, ub, grad, xresv, gradTmp;
};











