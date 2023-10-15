#pragma once


#ifndef vec
#define vec std::vector
#endif


struct Deltas
{
  int dim, size;
  vec<char> v;
  Deltas(int dim) { reset(dim); }
  Deltas() {}
  void reset(int dim)
  {
    this->dim = dim;
    size = 1;
    for (int i = 0; i < dim; ++i) size *= 3; // -1, 0, 1
    v.resize(size * dim);
    for (int i = 0; i < dim; ++i) v[i] = -1;
    char *prior = &v[0];
    for (int k = 1; k < size; ++k)
    {
      char *current = prior + dim;
      std::copy(prior, prior + dim, current);
      current[0] += 1;
      for (int i = 0, iend = dim - 1; i < iend; ++i)
      {
        if (current[i] <= 1) break;
        current[i] = -1;
        current[i + 1] += 1;
      }
      prior = current;
    }
  }
  
  
  char *operator()(int k)
  {
    return &v[0] + k * dim;
  }
  
  
  int Nvec() { return size; }
  int vecDim() { return dim; }
  
  
};


template<int dim>
struct GridID 
{ 
  int x[dim]; 
  int & operator[](int i) { return x[i]; }
};


template<int dim, typename itype>
struct ArrayHash
{
  std::size_t operator()(std::array<dim, itype> &x)
  {
    
  }
};


template <int dim, typename gridIDint> // Problem dimension.
struct GridSearch
{
  int maxEval;
  Deltas deltas;
  std::unordered_set<std::array<gridIDint, dim> > S;
  
  
  GridSearch(int maxEval): maxEval(maxEval)
  {
    deltas.reset(dim);
    S.reserve(maxEval);
  }
  
  
  template<typename Fun, typename LR> // Functor and Learning rate functor
  void operator()(double *x, Fun &f, LR &lr)
  {
    double learingRates[dim];
    std::array<gridIDint, dim> id, idnew;
    std::fill(id.begin(), id.end(), 0);
    double coor[dim];
    double fvalMin = 1e300;
    
    
    for (int i = 0; i < maxEval; ++i)
    {
      lr(x, learingRates);
      int kmin = -1;
      for (int k = 0, kend = deltas.Nvec(); k < kend; ++k)
      {
        auto *d = deltas(k);
        for (int j = 0; j < dim; ++j) idnew[j] = id[j] + d[j];
        if (S.find(idnew) == S.end()) continue;
        
        
        for (int j = 0; j < dim; ++j) 
          coor[j] = x[j] + learingRates[j] * d[j];
        
        
        double fval = f(coor);
        if (fval < fvalMin) kmin = k;
        S.insert(idnew);
      }
      if (kmin == -1) break;
      
      
      auto *d = deltas(kmin);
      for (int j = 0; j < dim; ++j)
      {
        id[j] += d[j];
        x[j] += learingRates[j] * d[j];
      }
    }
  }
  
};















// // Simple gradient descent.
// template<typename ing, typename num>
// struct GD
// {
//   ing maxit;
//   num lr; // Learning rate
//   GD(ing maxit, num lr): maxit(maxit), lr(lr) {}
//   
//   
//   template <typename Fun>
//   void operator()(num *x, ing dim, Fun f)
//   {
//     num priorVal = 1e38;
//     for (ing i = 0; i < maxit; ++i)
//     {
//       num val = f(x, dim);
//       num *grad = f.grad();
//       
//     }
//     
//     
//   }
//   
// };


#ifdef vec
#undef vec
#endif








