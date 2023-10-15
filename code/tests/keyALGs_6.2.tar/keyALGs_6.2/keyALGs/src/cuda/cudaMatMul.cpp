// #include <iostream>
// #include <fstream>
#include <stdio.h>
// #include <random>
#define gpuCompute


template<typename ing>
#ifdef gpuCompute
__device__
#endif
void jobIDtoParameter(ing id, ing NrowBlockShift, ing rowBsize,
                      ing &whichYcol, ing &rowStart)
{
  whichYcol = id >> NrowBlockShift;
  rowStart = (id & ((1 << NrowBlockShift) - 1)) * rowBsize;
}


template<typename ing, typename num>
#ifdef gpuCompute
__device__
#endif
void core(num *X, num *Y, num *Z,
          ing xrow, ing xcol, ing whichYcol,
          ing rowStart, ing rowBsize)
{
  ing &yrow = xcol, &zrow = xrow;
  ing xrowEnd = rowStart + rowBsize >= xrow ? xrow : rowStart + rowBsize;
  num *y = Y + yrow * whichYcol, *z = Z + zrow * whichYcol;
  for(ing i = 0; i < xcol; ++i)
  {
    num *xi = X + i * xrow;
    for(ing j = rowStart; j < xrowEnd; ++j) z[j] += xi[j] * y[i];
  }
}


#ifdef gpuCompute
template<typename ing, typename num>
__global__ void kernel(ing Njobs, num *X, num *Y, num *Z,
                       ing xrow, ing xcol, ing rowBsize, ing NrowBlockShift)
{
  ing jobID = blockIdx.x * blockDim.x + threadIdx.x;
  if(jobID >= Njobs) return;
  ing whichYcol, rowStart;
  jobIDtoParameter<ing>(jobID, NrowBlockShift, rowBsize, whichYcol, rowStart);
  core(X, Y, Z, xrow, xcol, whichYcol, rowStart, rowBsize);
}
#endif


#ifndef gpuCompute
template<typename ing, typename num>
void kernel(ing Njobs, ing jobID, num *X, num *Y, num *Z,
            ing xrow, ing xcol, ing rowBsize, ing NrowBlockShift)
{
  if(jobID >= Njobs) return;
  ing whichYcol, rowStart;
  jobIDtoParameter<ing>(jobID, NrowBlockShift, rowBsize, whichYcol, rowStart);
  core(X, Y, Z, xrow, xcol, whichYcol, rowStart, rowBsize);
}
#endif


template<typename T>
T Max(T a, T b) { return a > b ? a : b; }
template<typename ing, typename num>
void jobParameter(ing xrow, ing xcol, ing ycol, ing &rowBsize,
                  ing &NrowBlock, ing &NrowBlockShift, ing &NtotalJobs)
{
  NrowBlock = (xrow + rowBsize - 1) / rowBsize;
  // Make the number of row blocks power of 2:
  ing i = 1;
  NrowBlockShift = 0;
  while(i <= NrowBlock) { i *= 2; ++NrowBlockShift; }
  i = Max<ing> (1, i / 2);
  NrowBlockShift = Max<ing> (0, NrowBlockShift - 1);
  NrowBlock = i;
  rowBsize = (xrow + NrowBlock - 1) / NrowBlock;
  NtotalJobs = NrowBlock * ycol;
}


template<typename ing, typename num>
void matmul(num *X, num *Y, num *Z, ing xrow, ing xcol, ing ycol, ing rowBsize)
{
  ing NrowBlock, NtotalJobs, NrowBlockShift;
  jobParameter<ing, num>(xrow, xcol, ycol, rowBsize, NrowBlock,
                         NrowBlockShift, NtotalJobs);


#ifdef gpuCompute
  num *cudaX, *cudaY, *cudaZ;
  cudaMalloc(&cudaX, sizeof(num) * xrow * xcol);
  cudaMemcpy(cudaX, X, sizeof(num) * xrow * xcol, cudaMemcpyHostToDevice);
  cudaMalloc(&cudaY, sizeof(num) * xcol * ycol);
  cudaMemcpy(cudaY, Y, sizeof(num) * xcol * ycol, cudaMemcpyHostToDevice);
  cudaMalloc(&cudaZ, sizeof(num) * xrow * ycol);
  cudaMemset(cudaZ, 0, sizeof(num) * xrow * ycol);
  kernel<ing, num><<<(NtotalJobs + 256 - 1) / 256, 256>>> (
      NtotalJobs, cudaX, cudaY, cudaZ, xrow, xcol, rowBsize, NrowBlockShift);
  cudaMemcpy(Z, cudaZ, sizeof(num) * xrow * ycol, cudaMemcpyDeviceToHost);
  cudaFree(cudaX);
  cudaFree(cudaY);
  cudaFree(cudaZ);
#else
  for(ing i = 0; i < NtotalJobs; ++i)
    kernel<ing, num>(NtotalJobs, i, X, Y, Z, xrow, xcol,
                     rowBsize, NrowBlockShift);
  // for(ing i = 0, iend = 20; i < iend; ++i) Z[i] = 1;
#endif
  // printf("NtotalJobs = %d\n", int(NtotalJobs));
   // << "NtotalJobs = " << NtotalJobs << "\n";
}




#ifdef __cplusplus
extern "C" {
#endif

  void __declspec(dllexport)
  matmuldoubleCharlie(
      double *X, double *Y, double *Z,
      int *xrow, int *xcol, int *ycol, int *rowBsize
  );

#ifdef __cplusplus
}
#endif




void matmuldoubleCharlie(
    double *X, double *Y, double *Z,
    int *xrow, int *xcol, int *ycol, int *rowBsize
)
{
  matmul<size_t, double>(X, Y, Z, *xrow, *xcol, *ycol, *rowBsize);
}





/*
 template<typename ing, typename num>
 void testMatmul(num *X, num *Y, num *Z,
 ing xrow, ing xcol, ing ycol,
 ing rowBsize)
 {
 ing NrowBlock, NtotalJobs;
 jobParameter<ing, num> (xrow, xcol, ycol, rowBsize, NrowBlock, NtotalJobs);
 for(ing id = 0; id < NtotalJobs; ++id)
 {
 ing whichYcol, rowStart;
 jobIDtoParameter<ing, num> (id, NrowBlock, rowBsize, whichYcol, rowStart);
 core(X, Y, Z, xrow, xcol, whichYcol, rowStart, rowBsize);
 }
 }


 // [[Rcpp::export]]
 NumericMatrix testMatMulRcpp(NumericMatrix X, NumericMatrix Y, int rowBsize = 3)
 {
 NumericMatrix Z(X.nrow(), Y.ncol());
 testMatmul<int, double> (&X[0], &Y[0], &Z[0],
 X.nrow(), X.ncol(), Y.ncol(), rowBsize);
 return Z;
 }
 */












