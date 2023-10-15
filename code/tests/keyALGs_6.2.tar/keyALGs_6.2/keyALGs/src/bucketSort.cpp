// #pragma once


#include <Rcpp.h>
#include <chrono>
using namespace Rcpp;
#define vec std::vector


// Style is for C code.
char *ptr2content;
char **ptr2contentRef = &ptr2content;
void assignPtr2content(void *p) { *ptr2contentRef = (char*)p; }
int contentEleByteSize;
int *contentEleByteSizeRef = &contentEleByteSize;
void assignContentEleByteSize(int integer) { *contentEleByteSizeRef = integer; }
int sharedCmp(const void *a, const void *b)
{
  // a and b are essentially pointers to integers.
  return memcmp(ptr2content + *((int*)a) * contentEleByteSize,
                ptr2content + *((int*)b) * contentEleByteSize,
                contentEleByteSize);
}


void qsortInd(int *id, int Nele, void *X, int XeleByte)
{
  assignPtr2content(X);
  assignContentEleByteSize(XeleByte);
  qsort(id, Nele, sizeof(int), sharedCmp);
}


// Will only allow scan 1 or 2 bytes.
/*
 vec<int> sampleSort(char *X, size_t XeleByte, int startByte,
 int endByte, int *id, int Nelement)
 {
 endByte = endByte < (int)XeleByte ? endByte : XeleByte;
 int keySize = endByte - startByte;
 if(keySize <= 0) return vec<int> (0);
 keySize = keySize > 2 ? 2 : keySize;


 int containerSize = keySize == 2 ? 256 * 256 : 256;


 vec<vec<int> > C(containerSize);
 if(keySize == 1)
 {
 for(int i = 0; i < Nelement; ++i)
 {
 unsigned char k = *((unsigned char*)(X + id[i] * XeleByte + startByte));
 C[k].push_back(id[i]);
 }
 }
 else
 {
 for(int i = 0; i < Nelement; ++i)
 {
 unsigned char *k = (unsigned char*)(X + id[i] * XeleByte + startByte);
 unsigned short m = 256;
 C[k[0] * m + k[1]].push_back(id[i]);
 }
 }


 vec<int> rst(Nelement); rst.resize(0);


 for(int i = 0; i < containerSize; ++i)
 {
 int CiSize = C[i].size();
 if(CiSize == 0) continue;
 if(CiSize == 1) { rst.push_back(C[i][0]); continue; }


 int *copyStart = &*rst.end();
 rst.resize(rst.size() + CiSize);
 if(CiSize <= 500)
 {
 qsortInd(&C[i][0], CiSize, X, XeleByte);
 memcpy(copyStart, &C[i][0], CiSize * sizeof(int));
 }
 else
 {
 vec<int> subrst = sampleSort(
 X, XeleByte, endByte, endByte + 1, &C[i][0], CiSize);


 if(subrst.size() != 0)
 memcpy(copyStart, &subrst[0], CiSize * sizeof(int));
 else memcpy(copyStart, &C[i][0], CiSize * sizeof(int));
 }
 }


 return rst;
 }
 */


typedef struct { int size, capacity, *val; } Seg;


// A recursive function.
int bucketSortCore(char *X, size_t XeleByte, int startByte,
                   int endByte, int *id, int Nelement, int *rst)
{
  endByte = endByte < (int)XeleByte ? endByte : XeleByte;
  int keySize = endByte - startByte;
  if(keySize <= 0) return 0;
  keySize = keySize > 2 ? 2 : keySize;
  endByte = startByte + keySize;


  int Nbuckets = keySize == 2 ? 256 * 256 : 256;
  Seg bucket[Nbuckets];
  memset(bucket, 0, sizeof(Seg) * Nbuckets);


  char *y = X + startByte;
  if(keySize == 1)
  {
    for(int i = 0; i < Nelement; ++i)
    {
      unsigned char k = *((unsigned char*)(y + id[i] * XeleByte));
      ++bucket[k].capacity;
    }
  }
  else
  {
    for(int i = 0; i < Nelement; ++i)
    {
      unsigned char *k = (unsigned char*)(y + id[i] * XeleByte);
      unsigned short m = 256;
      ++bucket[k[0] * m + k[1]].capacity;
    }
  }


  bucket[0].val = rst;
  for(int i = 1; i < Nbuckets; ++i)
    bucket[i].val = bucket[i - 1].val + bucket[i - 1].capacity;


  if(keySize == 1)
  {
    for(int i = 0; i < Nelement; ++i)
    {
      unsigned char k = *((unsigned char*)(y + id[i] * XeleByte));
      Seg *s = &bucket[k];
      s->val[s->size] = id[i];
      ++s->size;
    }
  }
  else
  {
    for(int i = 0; i < Nelement; ++i)
    {
      unsigned char *k = (unsigned char*)(y + id[i] * XeleByte);
      unsigned short m = 256;
      Seg *s = &bucket[k[0] * m + k[1]];
      s->val[s->size] = id[i];
      ++s->size;
    }
  }


  int didsort = 1;
  for(int i = 0; i < Nbuckets; ++i)
  {
    if( bucket[i].size <= 1 ) continue;
    if ( bucket[i].size <= 500 )
    {
      qsortInd ( bucket[i].val, bucket[i].size, X, XeleByte );
    }
    else
    {
      int *subrst = id + (bucket[i].val - rst);
      didsort = bucketSortCore(X, XeleByte, endByte, endByte + 1,
                               bucket[i].val, bucket[i].size, subrst);
      if(didsort) memcpy(bucket[i].val, subrst, sizeof(int) * bucket[i].size);
    }
  }


  return didsort;
}


// Function returns a pointer to the result living in buffer B.
// B is an existing buffer of size >= 2 * Nelement in terms of int.
int *bucketSort(char *X, size_t XeleByte, int Nelement, int *B)
{
  int *id = B, *rst = id + Nelement;
  for(int i = 0; i < Nelement; ++i) id[i] = i;
  bucketSortCore(X, XeleByte, 0, 2, id, Nelement, rst);
  return rst;
}


// [[Rcpp::export]]
void testSampleSort(IntegerVector x, int eleByteSize = 4)
{
  int N = x.size() * sizeof(int) / eleByteSize;
  char *y = (char*)&x[0];
  vec<int> ind(N);
  std::iota(ind.begin(), ind.end(), 0);


  IntegerVector B(N * 2);
  auto now = std::chrono::steady_clock::now();
  int *bucketSortRst = bucketSort(y, eleByteSize, N, &B[0]);
  auto timecost = std::chrono::duration_cast<std::chrono::microseconds> (
    std::chrono::steady_clock::now() - now).count();
  Rcout << "sample sort time = " << timecost << "\n";


  now = std::chrono::steady_clock::now();
  qsortInd(&ind[0], N, y, eleByteSize);
  timecost = std::chrono::duration_cast<std::chrono::microseconds> (
    std::chrono::steady_clock::now() - now).count();
  Rcout << "quick sort time = " << timecost << "\n";


  for(int i = 1; i < N; ++i)
  {
    char *previous = y + bucketSortRst[i - 1] * eleByteSize;
    char *current = y + bucketSortRst[i] * eleByteSize;
    int tmp = memcmp(previous, current, eleByteSize);
    if(tmp > 0) { Rcout << "Sample sort not sorted"; break; }
  }


  for(int i = 1, iend = ind.size(); i < iend; ++i)
  {
    char *previous = y + ind[i - 1] * eleByteSize;
    char *current = y + ind[i] * eleByteSize;
    int tmp = memcmp(previous, current, eleByteSize);
    if(tmp > 0) { Rcout << "Quick sort not sorted"; break; }
  }
}


// [[Rcpp::export]]
void testSampleSortOnKeys(IntegerMatrix key)
{
  int eleByteSize = key.nrow();
  int N = key.ncol();
  IntegerVector v((N * eleByteSize + 3) / 4);
  unsigned char *x = (unsigned char*)(&v[0]);
  for(int i = 0, iend = N * eleByteSize; i < iend; ++i) x[i] = key[i];
  char *y = (char*)x;


  vec<int> ind(N);
  std::iota(ind.begin(), ind.end(), 0);
  IntegerVector B(2 * N);
  auto now = std::chrono::steady_clock::now();
  // sampleSort(y, eleByteSize, 0, 2, &ind[0], N, &sampleSortRst[0]);
  int *bucketSortRst = bucketSort(y, eleByteSize, N, &B[0]);
  auto timecost = std::chrono::duration_cast<std::chrono::microseconds> (
    std::chrono::steady_clock::now() - now).count();
  Rcout << "sample sort time = " << timecost << "\n";



  now = std::chrono::steady_clock::now();
  qsortInd(&ind[0], N, y, eleByteSize);
  timecost = std::chrono::duration_cast<std::chrono::microseconds> (
    std::chrono::steady_clock::now() - now).count();
  Rcout << "quick sort time = " << timecost << "\n";


  for(int i = 1; i < N; ++i)
  {
    char *previous = y + bucketSortRst[i - 1] * eleByteSize;
    char *current = y + bucketSortRst[i] * eleByteSize;
    int tmp = memcmp(previous, current, eleByteSize);
    if(tmp > 0) { Rcout << "Sample sort not sorted"; break; }
  }


  for(int i = 1, iend = ind.size(); i < iend; ++i)
  {
    char *previous = y + ind[i - 1] * eleByteSize;
    char *current = y + ind[i] * eleByteSize;
    int tmp = memcmp(previous, current, eleByteSize);
    if(tmp > 0) { Rcout << "Quick sort not sorted"; break; }
  }
}


inline int extractIntFromByte(char *x, int startByte, int endByte)
{
  unsigned char *y = (unsigned char*)x;
  int rst = 0, shift = 0;
  for(int i = endByte - 1; i >= startByte; --i, shift += 8)
    rst += ((int)(y[i]) << shift);
  return rst;
}


// For every element X[?], find the first element in Y such that
// it is no less than X[?].
// Xid should be in [0, NXele), Yid should be in [0, NYele).


/*
 vec<std::pair<int, int> > lowerBounds(
 char *X, char *Y, int eleByte, int startByte, int endByte,
 int *Xid, int *Yid, int NXele, int NYele, int &counter)
 {
 endByte = endByte < (int)eleByte ? endByte : eleByte;
 int keySize = endByte - startByte;
 if(keySize <= 0) return vec<std::pair<int, int> > (0);
 keySize = keySize > 2 ? 2 : keySize;
 endByte = startByte + keySize;


 int Nbuckets = keySize == 2 ? 256 * 256 : 256;


 vec<vec<int> > Xbucket(Nbuckets), Ybucket(Nbuckets);
 for(int i = 0; i < NXele; ++i)
 {
 int k = Xid[i];
 int l = extractIntFromByte(X + k * eleByte, startByte, endByte);
 Xbucket[l].push_back(k);
 }


 for(int i = 0; i < NYele; ++i)
 {
 int k = Yid[i];
 int l = extractIntFromByte(Y + k * eleByte, startByte, endByte);
 Ybucket[l].push_back(k);
 }


 // Rcout << "Xbucket = \n";
 // for(int i = 0, iend = Xbucket.size(); i < iend; ++i)
 // {
 //   if(Xbucket[i].size() != 0)
 //   {
 //     Rcout << "i = " << i << "\n";
 //     for(int k = 0, kend = Xbucket[i].size(); k < kend; ++k)
 //       Rcout << Xbucket[i][k] << ", ";
 //     Rcout << "\n";
 //   }
 // }


 // Rcout << "Ybucket = \n";
 // for(int i = 0, iend = Ybucket.size(); i < iend; ++i)
 // {
 //   if(Ybucket[i].size() != 0)
 //   {
 //     Rcout << "i = " << i << "\n";
 //     for(int k = 0, kend = Ybucket[i].size(); k < kend; ++k)
 //       Rcout << Ybucket[i][k] << ", ";
 //     Rcout << "\n";
 //   }
 // }
 // Rcout << "\n\n";




 vec<std::pair<int, int> > rst; rst.reserve(NXele);
 for(int i = 0, iend = Ybucket.size(); i < iend; ++i)
 {
 if(Ybucket[i].size() == 0)
 {
 for(int k = 0, kend = Xbucket[i].size(); k < kend; ++k)
 {
 rst.push_back(std::pair<int, int> (Xbucket[i][k], counter));
 }
 }
 else
 {
 if(Xbucket[i].size() == 0)
 {
 // counter += Ybucket[i].size();
 counter += 1;
 }
 else
 {
 // char *X, char *Y, int eleByte, int startByte, int endByte,
 // int *Xind, int *Yind, int NXele, int NYele
 if(eleByte == endByte) // We ARE comparing the last bytes of X and Y.
 {
 for(int k = 0, kend = Xbucket[i].size(); k < kend; ++k)
 rst.push_back(std::pair<int, int> (Xbucket[i][k], counter));
 counter += 1;
 }
 else
 {
 vec<std::pair<int, int> > subrst = lowerBounds(
 X, Y, eleByte, endByte, endByte + 1, &Xbucket[i][0], &Ybucket[i][0], Xbucket[i].size(),
 Ybucket[i].size(), counter);
 rst.resize(rst.size() + subrst.size());
 std::copy(subrst.begin(), subrst.end(), rst.begin() + rst.size() - subrst.size());
 // counter += Ybucket[i].size();
 }
 }
 }
 }


 return rst;
 }
 */


// typedef struct { int sizeX, capacityX, sizeY, capacityY, *x, *y; } Sxyeg;


typedef struct { int first, second; } Pair;


void bucketlb(
    char *X, char *Y, int eleByte, int startByte, int endByte,
    int *Xid, int *Yid, int NXele, int NYele, int *counter,
    int *Xbuffer, int *Ybuffer, Pair *rst)
{
  endByte = endByte < (int)eleByte ? endByte : eleByte;
  int keySize = endByte - startByte;
  // if(keySize <= 0) return vec<std::pair<int, int> > (0);
  if(keySize <= 0) return;
  keySize = keySize > 2 ? 2 : keySize;
  endByte = startByte + keySize;
  int Nbuckets = keySize == 2 ? 256 * 256 : 256;


  Seg B[Nbuckets * 2];
  memset(B, 0, sizeof(Seg) * Nbuckets * 2);
  Seg *Xbucket = B, *Ybucket = Xbucket + Nbuckets;


  for(int i = 0; i < NXele; ++i)
  {
    int l = extractIntFromByte(X + Xid[i] * eleByte, startByte, endByte);
    ++Xbucket[l].capacity;
  }
  Xbucket[0].val = Xbuffer;
  for(int i = 1; i < Nbuckets; ++i)
    Xbucket[i].val = Xbucket[i - 1].val + Xbucket[i - 1].capacity;
  for(int i = 0; i < NXele; ++i)
  {
    int l = extractIntFromByte(X + Xid[i] * eleByte, startByte, endByte);
    Seg *s = &Xbucket[l];
    s->val[s->size] = Xid[i];
    ++s->size;
  }


  for(int i = 0; i < NYele; ++i)
  {
    int l = extractIntFromByte(Y + Yid[i] * eleByte, startByte, endByte);
    ++Ybucket[l].capacity;
  }
  Ybucket[0].val = Ybuffer;
  for(int i = 1; i < Nbuckets; ++i)
    Ybucket[i].val = Ybucket[i - 1].val + Ybucket[i - 1].capacity;
  for(int i = 0; i < NYele; ++i)
  {
    int l = extractIntFromByte(Y + Yid[i] * eleByte, startByte, endByte);
    Seg *s = &Ybucket[l];
    s->val[s->size] = Yid[i];
    ++s->size;
  }


  // Rcout << "Xbucket = \n";
  // for(int i = 0, iend = Nbuckets; i < iend; ++i)
  // {
  //   if(Xbucket[i].size != 0)
  //   {
  //     Rcout << "i = " << i << "\n";
  //     for(int k = 0, kend = Xbucket[i].size; k < kend; ++k)
  //       Rcout << Xbucket[i].val[k] << ", ";
  //     Rcout << "\n";
  //   }
  // }


  // Rcout << "Ybucket = \n";
  // for(int i = 0, iend = Nbuckets; i < iend; ++i)
  // {
  //   if(Ybucket[i].size != 0)
  //   {
  //     Rcout << "i = " << i << "\n";
  //     for(int k = 0, kend = Ybucket[i].size; k < kend; ++k)
  //       Rcout << Ybucket[i].val[k] << ", ";
  //     Rcout << "\n";
  //   }
  // }
  // Rcout << "\n\n";


  // Xid and Yid are now useless.
  int rstI = 0;
  for(int i = 0; i < Nbuckets; ++i)
  {
    if(Ybucket[i].size == 0)
    {
      for(int k = 0, kend = Xbucket[i].size; k < kend; ++k)
      {
        rst[rstI].first = Xbucket[i].val[k];
        rst[rstI].second = *counter;
        // Rcout << "first, second = " << rst[k].first
        // << ", " << rst[k].second << "\n";
        ++rstI;
      }
      // Rcout << "\n";
    }
    else
    {

      // Rcout << "counter before = " << counter << ", ";


      if(Xbucket[i].size == 0) *counter += 1;
      else
      {
        if(eleByte == endByte) // We ARE comparing the last bytes of X and Y.
        {
          for(int k = 0, kend = Xbucket[i].size; k < kend; ++k)
          {
            rst[rstI].first = Xbucket[i].val[k];
            rst[rstI].second = *counter;
            // Rcout << "first, second = " << rst[k].first
            // << ", " << rst[k].second << "\n";
            ++rstI;
          }
          // Rcout << "\n";
          *counter += 1;
        }
        else
        {
          // lowerBounds(
          //   char *X, char *Y, int eleByte, int startByte, int endByte,
          //   int *Xid, int *Yid, int NXele, int NYele, int &counter,
          //   int *Xbuffer, int *Ybuffer, Pair *rst)
          bucketlb(
            X, Y, eleByte, endByte, endByte + 1,
            Xbucket[i].val, Ybucket[i].val, Xbucket[i].size, Ybucket[i].size, counter,
            Xid + (Xbucket[i].val - Xbuffer), Yid + (Ybucket[i].val - Ybuffer),
            rst + (Xbucket[i].val - Xbuffer)  );
          rstI += Xbucket[i].size;
        }
      }


      // Rcout << "counter after = " << counter << "\n";

    }


  }


}


// Function returns the pointer to the result that is stored in B. The result
// will be of size NXele.
// B is an existing buffer of size >= 4 * NXele + 2 * NYele, in terms of int.
int *bucketingLB(char *X, char *Y, int eleByte, int NXele, int NYele, int *B)
{
  int *Xid = B, *Yid = B + NXele;
  int *Xbuffer = Yid + NYele, *Ybuffer = Xbuffer + NXele;
  Pair *p = (Pair*)(Ybuffer + NYele);
  for(int i = 0; i < NXele; ++i) Xid[i] = i;
  for(int i = 0; i < NYele; ++i) Yid[i] = i;
  int *rst = (int*)(p) - NXele;
  int counter = 0;
  bucketlb(X, Y, eleByte, 0, 1, Xid, Yid, NXele, NYele,
           &counter, Xbuffer, Ybuffer, p);
  for(int i = 0; i < NXele; ++i) rst[p[i].first] = p[i].second;
  return rst;
}


vec<char> intArray2charArray(int *x, int size)
{
  vec<char> rst; rst.reserve(size * sizeof(int));
  for(int i = 0; i < size; ++i)
  {
    char *p = (char*)(x + i);
    for(int k = sizeof(int) - 1; k >= 0; --k) rst.push_back(p[k]);
  }
  return rst;
}


int lower_bound(char *arr, int arrSize, char *keyval, int eleByte)
{
  int mid, low = 0, high = arrSize, cmpRst;
  while (low < high)
  {
    mid = low + (high - low) / 2;
    cmpRst = memcmp(keyval, arr + mid * eleByte, eleByte);
    if (cmpRst <= 0) high = mid;
    else low = mid + 1;
  }
  cmpRst = memcmp(arr + low * eleByte, keyval, eleByte);
  if(low < arrSize && cmpRst < 0) ++low;
  return low;
}


// [[Rcpp::export]]
IntegerVector testLowerBounds(IntegerVector x, IntegerVector y, int eleByte = 4)
{
  vec<char> Cx = intArray2charArray(&x[0], x.size());
  vec<char> Cy = intArray2charArray(&y[0], y.size());
  int NXele = Cx.size() / eleByte, NYele = Cy.size() / eleByte;
  vec<int> buffer(NXele * 4 + 2 * NYele);
  int *rst = bucketingLB(&Cx[0], &Cy[0], eleByte, NXele, NYele, &buffer[0]);
  return IntegerVector(rst, rst + NXele) + 1;
}


// [[Rcpp::export]]
void testLowerBoundsOnKeys(IntegerMatrix key, IntegerVector splitterInd)
{
  int eleByteSize = key.nrow();
  int N = key.ncol();
  IntegerVector v((N * eleByteSize + 3) / 4);
  unsigned char *x = (unsigned char*)(&v[0]);
  for(int i = 0, iend = N * eleByteSize; i < iend; ++i) x[i] = key[i];
  char *y = (char*)x;


  int Nsplitter = splitterInd.size();
  vec<char> splitter(Nsplitter * eleByteSize);
  for(int i = 0; i < Nsplitter; ++i)
  {
    memcpy(&splitter[0] + i * eleByteSize, splitterInd[i] * eleByteSize + &key[0], eleByteSize);
  }


  IntegerVector lbRst(N);
  auto now = std::chrono::steady_clock::now();
  for(int i = 0; i < N; ++i)
    lbRst[i] = lower_bound(y, N, &splitter[0], eleByteSize);
  auto timecost = std::chrono::duration_cast<std::chrono::microseconds> (
    std::chrono::steady_clock::now() - now).count();
  Rcout << "lower bound time = " << timecost << "\n";


  IntegerVector B(4 * N + 2 * Nsplitter);
  now = std::chrono::steady_clock::now();
  bucketingLB(y, &splitter[0], eleByteSize, N, Nsplitter, &B[0]);
  timecost = std::chrono::duration_cast<std::chrono::microseconds> (
    std::chrono::steady_clock::now() - now).count();
  Rcout << "bucket lower bound time = " << timecost << "\n";
}






























