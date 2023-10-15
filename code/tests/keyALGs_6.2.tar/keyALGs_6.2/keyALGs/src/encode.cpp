// #pragma once
// ============================================================================
// # include <setjmp.h>
// jmp_buf env;
// longjmp(env, 1); // This line is the breaking point.
// if(setjmp(env)) return List::create(); // Put this line at the beginning of
// the export function.
// ============================================================================
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
// # define ARMA_DONT_PRINT_ERRORS
// // [[Rcpp::depends(RcppArmadillo)]]
// # include <RcppArmadillo.h>
# include <Rcpp.h>
// # include <RcppParallel.h>
// # include <chrono>
# include <fstream>
# include <iostream>
# include <bitset>
// # include <random>
// # include "h/dnyTasking.hpp"
# include "pcg/pcg_random.hpp"
# define vec std::vector
# define RNG pcg32
using namespace Rcpp;


void threeBytes2fourChar(std::string &lookup,
            char *threebytes, int threebytesN,
            char *fourChars)
{
  int x = 0;
  std::copy(threebytes, threebytes + threebytesN, (char*)(&x));


  fourChars[0] = lookup[x & 0b111111];
  fourChars[1] = lookup[(x & 0b111111000000) >> 6];
  fourChars[2] = lookup[(x & 0b111111000000000000) >> 12];
  fourChars[3] = lookup[(x & 0b111111000000000000000000) >> 18];

}


void fourChar2threeBytes(std::unordered_map<char, unsigned char> &lookup,
                         char *threebytes, int threebytesN, char *fourChars)
{
  int x = int(lookup[fourChars[0]]) | (int(lookup[fourChars[1]]) << 6) |
    (int(lookup[fourChars[2]]) << 12) | (int(lookup[fourChars[3]]) << 18);


  std::copy((char*)(&x), (char*)(&x) + threebytesN, threebytes);
}




// [[Rcpp::export]]
void encode(std::string inputPath, std::string filetype, std::string outputTxtName)
{
  // Every 3 bytes are encoded by 4 characters.
  // First 8 bytes store the size of x.
  // Next 256 bytes store fileName.
  std::ifstream input(inputPath, std::ios::in | std::ios::binary);
  vec<char> x(std::istreambuf_iterator<char>(input), {});
  std::size_t xsize = x.size();
  std::string lookup = "WgP5D8w3IpRt ksBH2hLnzocGNYF-KC0iTQSOvl74ejmUJX19rEduyqfMaxAZ6bV";


  vec<char> y(8 + 256 + x.size());
  std::copy((char*)(&xsize), (char*)(&xsize) + 8, &y[0]);
  std::copy(filetype.begin(), filetype.end(), y.begin() + 8);
  std::copy(x.begin(), x.end(), y.begin() + 8 + 256);
  x.swap(y);
  vec<char>().swap(y);
  xsize = x.size();


  RNG rng(314159265);
  std::uniform_int_distribution<unsigned char> U(0, 255);
  for(auto & u: x) u ^= U(rng);


  std::size_t totalSize = (xsize / 3 + 1) * 4;
  std::string rst(totalSize, ' ');
  for(std::size_t j = 0, i = 0; j < xsize; i += 4, j += 3)
  {
    threeBytes2fourChar(lookup, &x[j], std::min<int64_t>(xsize - j, 3), &rst[i]);
  }


  std::ofstream out(outputTxtName);
  out << rst;
  out.close();
}




// [[Rcpp::export]]
void decode(std::string inputTxtPath, std::string outputNameNoFileType)
{

  std::ifstream input(inputTxtPath, std::ios::binary);
  vec<char> rst(std::istreambuf_iterator<char>(input), {});
  vec<char> theSize(12);
  std::copy(rst.begin(), rst.begin() + 12, theSize.begin());


  std::string previousLookup = "WgP5D8w3IpRt ksBH2hLnzocGNYF-KC0iTQSOvl74ejmUJX19rEduyqfMaxAZ6bV";
  std::unordered_map<char, unsigned char> lookup(64);
  for(int i = 0; i < 64; ++i) lookup[previousLookup[i]] = i;


  std::size_t xsize = 0;
  char *tmp0 = (char*)(&xsize);
  for(std::size_t j = 0, i = 0; j < 8; i += 4, j += 3)
  {
    fourChar2threeBytes(lookup, &tmp0[j], std::min<int64_t>(8 - j, 3), &theSize[i]);
  }


  RNG rng(314159265);
  std::uniform_int_distribution<unsigned char> U(0, 255);
  unsigned char *tmp1 = (unsigned char*)&xsize;
  for(int i = 0; i < 8; ++i) tmp1[i] ^= U(rng);


  std::string x(xsize + 256 + 8, ' ');
  xsize = x.size();


  for(std::size_t j = 0, i = 0; j < xsize; i += 4, j += 3)
  {
    fourChar2threeBytes(lookup, &x[j], std::min<int64_t>(xsize - j, 3), &rst[i]);
  }


  rng.seed(314159265);
  for(int i = 0, iend = x.size(); i < iend; ++i) x[i] ^= U(rng);


  std::string filetype;
  if(true)
  {
    int i = 8;
    for(int iend = 8 + 256; i < iend and x[i] != ' '; ++i);
    filetype.resize(i - 8);
    std::copy(&x[8], &x[i], &filetype[0]);
  }


  std::ofstream myFile(outputNameNoFileType + "." + filetype, std::ios::out | std::ios::binary);
  myFile.write(&x[8 + 256], xsize - (8 + 256));
}





































