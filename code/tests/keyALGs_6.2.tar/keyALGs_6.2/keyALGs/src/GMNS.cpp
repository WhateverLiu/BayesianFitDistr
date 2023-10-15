// #pragma once


// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
//#include <fast_mutex.h>
//#include <fstream>
#include <ctime>
using namespace RcppParallel;
using namespace Rcpp;


//tthread::mutex mutexm;


# define valtype double

struct G2{
valtype miu1, miu2, var1, var2, covar, w;
valtype*ptr;

valtype eval(valtype x1, valtype x2)
{
  valtype M=var1*var2-covar*covar;
  valtype tmp=w*std::exp((((x1-miu1)*(var2/M)+(x2-miu2)*(-covar/M))*(x1-miu1)+((x1-miu1)*
      (-covar/M)+(x2-miu2)*(var1/M))*(x2-miu2))/(-2))/(std::sqrt(M)*2*M_PI);
  if(std::isnan(tmp))return 0;
  return tmp;
}


valtype eigenRatioApp(){
return (var1+var2)/std::sqrt(var1*var2-covar*covar);
}


void evalV(valtype*x1, valtype*x2, int siz, valtype*rst, valtype eigenvalueRatioThreshold)
{
  valtype M=var1*var2-covar*covar;
  valtype eigenvalRatioApp=eigenRatioApp();
  valtype var2_M=var2/M, _covar_M=-covar/M, var1_M=var1/M, aboutPi=std::sqrt(M)*2*M_PI;
  for(valtype*x1end=x1+siz; x1!=x1end; ++x1,++x2,++rst)
  {
    if(eigenvalRatioApp>eigenvalueRatioThreshold)
    {
      *rst=0;
      continue;
    }
    valtype x1_miu1=*x1-miu1, x2_miu2=*x2-miu2;
    *rst=w*std::exp(((x1_miu1*var2_M+x2_miu2*_covar_M)*x1_miu1+(x1_miu1*
      _covar_M+x2_miu2*var1_M)*x2_miu2)/(-2))/aboutPi;
    if(!std::isfinite(*rst))*rst=0;
  }
}



valtype evalNoWeight(valtype x1, valtype x2)
{
  valtype M=var1*var2-covar*covar;
  valtype tmp=std::exp((((x1-miu1)*(var2/M)+(x2-miu2)*(-covar/M))*(x1-miu1)+((x1-miu1)*
      (-covar/M)+(x2-miu2)*(var1/M))*(x2-miu2))/(-2))/(std::sqrt(M)*2*M_PI);
  if(std::isnan(tmp))return 0;
  return tmp;
}




void formLocovar(valtype x, valtype y, G2*embeddedKernelV,
                 std::vector<unsigned>&ind, std::vector<valtype>&w){

valtype E[2], E2[3];
E[0]=0;E[1]=0;E2[0]=0;E2[1]=0;E2[2]=0;

for(unsigned i=0,iend=ind.size();i!=iend;++i)
{
  E2[0]+=w[i]*(embeddedKernelV[ind[i]].var1+
    embeddedKernelV[ind[i]].miu1*embeddedKernelV[ind[i]].miu1);
  E2[1]+=w[i]*(embeddedKernelV[ind[i]].covar+
    embeddedKernelV[ind[i]].miu1*embeddedKernelV[ind[i]].miu2);
  E2[2]+=w[i]*(embeddedKernelV[ind[i]].var2+
    embeddedKernelV[ind[i]].miu2*embeddedKernelV[ind[i]].miu2);

  E[0]+=w[i]*embeddedKernelV[ind[i]].miu1;
  E[1]+=w[i]*embeddedKernelV[ind[i]].miu2;
}

var1=E2[0]-E[0]*E[0];
covar=E2[1]-E[0]*E[1];
var2=E2[2]-E[1]*E[1];

miu1=x;
miu2=y;
}






valtype evalMahalanobisD(valtype x1, valtype x2){
valtype M=var1*var2-covar*covar;
valtype tmp=((x1-miu1)*(var2/M)+(x2-miu2)*(-covar/M))*(x1-miu1)+((x1-miu1)*
    (-covar/M)+(x2-miu2)*(var1/M))*(x2-miu2);
if(std::isnan(tmp))return 0;
return tmp;
}


};









valtype absRelativeDiff(G2&base, G2&x){
valtype tmp, diff=0;

tmp=std::abs(2*x.w/(x.w+base.w)-1);
if(std::isnan(tmp))tmp=0;
//else if(std::isinf(tmp))tmp=std::abs(x.w-base.w);
diff+=tmp;

tmp=std::abs(2*x.miu1/(x.miu1+base.miu1)-1);
if(std::isnan(tmp))tmp=0;
//else if(std::isinf(tmp))tmp=std::abs(x.miu1-base.miu1);
diff+=tmp;

tmp=std::abs(2*x.miu2/(base.miu2+x.miu2)-1);
if(std::isnan(tmp))tmp=0;
//else if(std::isinf(tmp))tmp=std::abs(x.miu2-base.miu2);
diff+=tmp;

tmp=std::abs(2*x.covar/(base.covar+x.covar)-1);
if(std::isnan(tmp))tmp=0;
//else if(std::isinf(tmp))tmp=std::abs(x.covar-base.covar);
diff+=tmp;

tmp=std::abs(2*x.var1/(base.var1+x.var1)-1);
if(std::isnan(tmp))tmp=0;
//else if(std::isinf(tmp))tmp=std::abs(x.var1-base.var1);
diff+=tmp;

tmp=std::abs(2*x.var2/(base.var2+x.var2)-1);
if(std::isnan(tmp))tmp=0;
//else if(std::isinf(tmp))tmp=1;
diff+=tmp;
return diff;
}




valtype absDiff(G2&base, G2&x){
//valtype tmp, diff=0;

return std::abs(x.w-base.w)+std::abs(x.miu1-base.miu1)+
  std::abs(x.miu2-base.miu2)+std::abs(x.var1-base.var1)+
  std::abs(x.var2-base.var2)+std::abs(x.covar-base.covar);
}




valtype evalDiff(G2&base, G2&x, bool useRelative){
if(useRelative)return absRelativeDiff(base, x);
return absDiff(base, x);
}



// // I is the Ith CPU
// void processSomeoneQuit(volatile bool&someoneQuit, int I,
//                           int NofCore, int*iStart, int*iCurrent, int*iEnd,
//                           std::vector<int>&indexV, signed char*haulted, bool eraseCurrent=0){
//
// if(someoneQuit)
// {
//   //std::cout<<-(int)I;
//   if(haulted[I]==0)haulted[I]=1;
//
//   if(I!=NofCore-1)while(someoneQuit);
//   else
//   {
//     //std::cout<<-(int)I;
//     while(true)
//     {
//       //std::cout<<-(int)I;
//       bool allHaulted=1;
//       for(int k=0;k!=NofCore;++k)
//       {
//         if(haulted[k]==0)
//         {
//           allHaulted=0;
//           break;
//         }
//       }
//
//       if(allHaulted)break;
//     }
//         // erase the elements
//     if(!eraseCurrent)
//     {
//       for(int k=I;k>=0;--k)
//         indexV.erase(indexV.begin()+iStart[k], indexV.begin()+iCurrent[k]);
//     }
//     else
//     {
//       for(int k=I;k>=0;--k)
//         indexV.erase(indexV.begin()+iStart[k], indexV.begin()+iCurrent[k]+1);
//     }
//
//     // reset current i
//     int avgtask=indexV.size()/NofCore;
//     if(avgtask!=0)
//     {
//       for(int k=0;k!=NofCore;++k)
//       {
//         iStart[k]=k*avgtask;
//         iCurrent[k]=iStart[k];// iEnd and iStart point to the same container
//       }
//       iEnd[NofCore-1]=indexV.size();
//     }
//     else
//     {
//       // iStart is of size NofCore+1
//       // first NofCore+1-indexV.size() are 0
//       // the next will be 1, 2, 3...
//       int*tmp=iStart+NofCore+1-indexV.size();
//       std::fill(iStart, tmp, 0);
//       for(int k=1,kend=indexV.size();k<=kend;++k)tmp[k-1]=k;
//       std::copy(iStart, iStart+NofCore, iCurrent);
//     }
//
//     // set back someoneQuit
//     someoneQuit=0;
//   }
//   if(iStart[I]==iEnd[I])haulted[I]=-1;// means the thread should die
//   else haulted[I]=0;
// }
//
// }




// struct dynamicTasking{
// volatile bool someoneQuit;
// std::vector<unsigned>indexV;
// std::vector<char>container;
// unsigned*iStart, *iEnd, *iCurrent;
// unsigned NofCore;
// signed char*haulted;
//
// void reset(unsigned NofCPU, unsigned NofAtomTask){
//
// someoneQuit=0;
//
// NofCore=NofCPU;
//
// indexV.resize(NofAtomTask);
// for(unsigned i=0;i!=NofAtomTask;++i)indexV[i]=i;
//
// container.resize(((2*NofCore+1)*sizeof(unsigned)+
//   NofCore*sizeof(signed char))/sizeof(char));
//
// iStart=(unsigned*)&container[0];
// iEnd=iStart+1;
// iCurrent=iEnd+NofCore;
// haulted=(signed char*)(iCurrent+NofCore);
//
// unsigned avgtask=NofAtomTask/NofCore;
// for(unsigned i=0;i!=NofCore;++i)
// {
//   iStart[i]=i*avgtask;
//   iCurrent[i]=iStart[i];
//   haulted[i]=0;
// }
// iStart[NofCore]=NofAtomTask;
// }
//
// dynamicTasking(unsigned NofCPU, unsigned NofAtomTask){
// reset(NofCPU, NofAtomTask);
// }
//
// dynamicTasking(){}
//
// void monitorManage(unsigned cpuI, bool currentTaskInProgress){
// if(someoneQuit)
// {
//   //std::cout<<-(int)I;
//   if(haulted[cpuI]==0)haulted[cpuI]=1;
//
//   if(cpuI!=NofCore-1)while(someoneQuit);
//   else
//   {
//     //std::cout<<-(int)I;
//     while(true)
//     {
//       //std::cout<<-(int)I;
//       bool allHaulted=1;
//       for(int k=0,kend=NofCore;k!=kend;++k)
//       {
//         if(haulted[k]==0)
//         {
//           allHaulted=0;
//           break;
//         }
//       }
//
//       if(allHaulted)break;
//     }
//     // erase the elements
//     if(!currentTaskInProgress)
//     {
//       for(int k=cpuI;k>=0;--k)
//         indexV.erase(indexV.begin()+iStart[k], indexV.begin()+iCurrent[k]);
//     }
//     else
//     {
//       for(int k=cpuI;k>=0;--k)
//         indexV.erase(indexV.begin()+iStart[k], indexV.begin()+iCurrent[k]+1);
//     }
//
//     // reset current i
//     int avgtask=indexV.size()/NofCore;
//     if(avgtask!=0)
//     {
//       for(int k=0,kend=NofCore;k!=kend;++k)
//       {
//         iStart[k]=k*avgtask;
//         iCurrent[k]=iStart[k];// iEnd and iStart point to the same container
//       }
//       iEnd[NofCore-1]=indexV.size();
//     }
//     else
//     {
//       // iStart is of size NofCore+1
//       // first NofCore+1-indexV.size() are 0
//       // the next will be 1, 2, 3...
//       unsigned*tmp=iStart+NofCore+1-indexV.size();
//       std::fill(iStart, tmp, 0);
//       for(int k=1,kend=indexV.size();k<=kend;++k)tmp[k-1]=k;
//       std::copy(iStart, iStart+NofCore, iCurrent);
//     }
//     // set back someoneQuit
//     someoneQuit=0;
//   }
//   if(iStart[cpuI]==iEnd[cpuI])haulted[cpuI]=-1;// means the thread should die
//   else haulted[cpuI]=0;
// }
// }
//
// bool nextTaskID(unsigned cpuI, unsigned&taskID)
// {
//   if(iCurrent[cpuI]==iEnd[cpuI])someoneQuit=1;
//   monitorManage(cpuI, false);
//   if(haulted[cpuI]==-1)return 0;
//   taskID=indexV[iCurrent[cpuI]];
//   ++iCurrent[cpuI];
//   //std::cout<<taskID<<" ";
//   return 1;
// }
// };







struct dynamicTasking{
tbb::mutex m;
std::vector<unsigned>indexV;
dynamicTasking(){}
void reset(unsigned NofAtomTask){
indexV.resize(NofAtomTask);
for(unsigned i=0,iend=NofAtomTask;i!=iend;++i)
  indexV[i]=i;
}
dynamicTasking(unsigned NofAtomTask){
reset(NofAtomTask);
}
bool nextTaskID(unsigned&taskID)
{
  bool rst=1;
  m.lock();
  if(indexV.size()==0)rst=0;
  else
  {
    taskID=indexV.back();
    indexV.pop_back();
  }
  m.unlock();
  return rst;
}
};












struct firstForEvalDensities: public Worker{
G2*&gmix;
//int&gmixSize;
double*&LongBegin;
double*&LatBegin;
int&NofLocations;
double&ellipseAxisRatioThreshold;

// int&NofCore;
// volatile bool&someoneQuit;// at first this is 0
// std::vector<int>&indexV;// a vector intilized as 0:(length(object)-1)
// int*&iStart;
// //int*&iEnd;// an array of size NofCore
// int*&iCurrent;
// signed char*&haulted;

dynamicTasking&dT;


//Here is how to initialize class member of references
firstForEvalDensities(
G2*&gmix,
//int&gmixSize,
double*&LongBegin, double*&LatBegin, int&NofLocations,
double&ellipseAxisRatioThreshold,
// int&NofCore,
// volatile bool&someoneQuit,
// std::vector<int>&indexV,
// int*&iStart,
// //int*&iEnd,
// int*&iCurrent,
// signed char*&haulted
dynamicTasking&dT):

gmix(gmix),
//gmixSize(gmixSize),
LongBegin(LongBegin),
LatBegin(LatBegin), NofLocations(NofLocations),
ellipseAxisRatioThreshold(ellipseAxisRatioThreshold),
// NofCore(NofCore),
// someoneQuit(someoneQuit),
// indexV(indexV),
// iStart(iStart),
// //iEnd(iEnd),
// iCurrent(iCurrent),
// haulted(haulted)
dT(dT)
{}

//when writing operator() function, it must use std::size_t!! Anything else won't pass!!
void operator()(std::size_t st, std::size_t end){
for(std::size_t I=st;I!=end;++I)
{
  for(;;)
  {
    unsigned objI=0;
    //if(!dT.nextTaskID(I, objI))break;
    if(!dT.nextTaskID(objI))break;
    //std::cout<<-(int)objI;
    gmix[objI].evalV(LongBegin, LatBegin, NofLocations, gmix[objI].ptr, ellipseAxisRatioThreshold);
  }
}
}
};

















struct secondForEvalWeightM: public Worker{
G2*&gmix;
int&gmixSize;
// int&NofLocations;
// int&NofCore;
dynamicTasking&dT;

//Here is how to initialize class member of references
secondForEvalWeightM(
G2*&gmix, int&gmixSize,
//int&NofLocations, int&NofCore,
dynamicTasking&dT
):
gmix(gmix), gmixSize(gmixSize),
//NofLocations(NofLocations),NofCore(NofCore),
dT(dT)
{}

//when writing operator() function, it must use std::size_t!! Anything else won't pass!!
void operator()(std::size_t st, std::size_t end){
for(std::size_t I=st;I!=end;++I)
{
  for(;;)
  {
    unsigned objI=0;
    if(!dT.nextTaskID(objI))
    {
      break;
    }
    valtype tmpS=0;
    for(int i=0,iend=gmixSize;i!=iend;++i)
      tmpS+=gmix[i].ptr[objI];
    for(int i=0,iend=gmixSize;i!=iend;++i)
    {
      gmix[i].ptr[objI]/=tmpS;
      if(!std::isfinite(gmix[i].ptr[objI]))gmix[i].ptr[objI]=0;
    }
  }

}
}
};













struct thirdForEvalGaussians: public Worker{
//G2*&gmix;
//int&gmixSize;
std::vector<G2>&gmix;
int&NofLocations;
//int&NofCore;
//unsigned char*&gmixIndicator;
double*&valBegin;
double&sumLoss;
double&weightEPS;
double*&Long;
double*&Lat;
valtype*&diff;
bool&useRelativeDiff;
dynamicTasking&dT;

//Here is how to initialize class member of references
thirdForEvalGaussians(
std::vector<G2>&gmix,
int&NofLocations,
//int&NofCore,
//unsigned char*&gmixIndicator,
double*&valBegin, double&sumLoss,
double&weightEPS,
double*&Long,
double*&Lat,
valtype*&diff,
bool&useRelativeDiff,
dynamicTasking&dT
):
gmix(gmix),
NofLocations(NofLocations),
//NofCore(NofCore),
//gmixIndicator(gmixIndicator),
valBegin(valBegin),
sumLoss(sumLoss), weightEPS(weightEPS),
Long(Long),
Lat(Lat),
diff(diff),
useRelativeDiff(useRelativeDiff),
dT(dT)
{}

//when writing operator() function, it must use std::size_t!! Anything else won't pass!!
void operator()(std::size_t st, std::size_t end){
for(std::size_t I=st;I!=end;++I)
{
  //if(someoneQuit)continue;

  // int avgtask=gmix.size()/NofCore;
  // int istart=avgtask*I, iend=istart+avgtask;
  // if(I==NofCore-1)iend=gmix.size();

  //for(int i=iend-1;i>=istart;--i)
  for(;;)
  {
    unsigned i=0;
    //if(!dT.nextTaskID(I, i))break;
    if(!dT.nextTaskID(i))break;

    G2 gmixiTmp=gmix[i];

    double*lossi=valBegin;
    for(int j=0;j!=NofLocations;++j,++lossi)
      gmix[i].ptr[j]*=*lossi;

    valtype Nk=std::accumulate(gmix[i].ptr, gmix[i].ptr+NofLocations, 0.0);
    gmix[i].w=Nk/sumLoss;

    if(gmix[i].w<weightEPS)continue;

    gmix[i].miu1=0;
    gmix[i].miu2=0;
    gmix[i].var1=0;
    gmix[i].var2=0;
    gmix[i].covar=0;
    for(int j=0;j!=NofLocations;++j)
    {
      gmix[i].miu1+=Long[j]*gmix[i].ptr[j];
      gmix[i].miu2+=Lat[j]*gmix[i].ptr[j];
    }
    gmix[i].miu1/=Nk;
    gmix[i].miu2/=Nk;

    for(int j=0;j!=NofLocations;++j)
    {
      gmix[i].var1+=(Long[j]-gmix[i].miu1)*(Long[j]-gmix[i].miu1)*gmix[i].ptr[j];
      gmix[i].var2+=(Lat[j]-gmix[i].miu2)*(Lat[j]-gmix[i].miu2)*gmix[i].ptr[j];
      gmix[i].covar+=(Long[j]-gmix[i].miu1)*(Lat[j]-gmix[i].miu2)*gmix[i].ptr[j];
    }

    gmix[i].var1/=Nk;
    gmix[i].var2/=Nk;
    gmix[i].covar/=Nk;

    diff[I]+=evalDiff(gmixiTmp, gmix[i], useRelativeDiff);
  }
}
}
};







// [[Rcpp::export]]
List gm2dParallel(NumericVector Long, NumericVector Lat, NumericVector val, NumericVector weight,
              NumericVector miu1, NumericVector miu2, NumericVector var1, NumericVector var2,
              NumericVector covar, double weightEPS, double convergeEPS,
              int maxit, int showProgress, int NofCore, bool convergeRelativeDiff=1,
              double ellipseAxisRatioThreshold=14.2){

int intialG=miu1.size();
//std::cout<<"intialG="<<intialG<<std::endl;
double sumLoss=std::accumulate(val.begin(),val.end(),0.0);
int NofLocations=Long.size();

std::vector<G2>gmix(intialG);

//std::vector<unsigned char>gmixIndicator(gmix.size());

//std::vector<std::vector<valtype> >wM(intialG, std::vector<valtype>(NofLocations));

for(int i=0,iend=gmix.size();i!=iend;++i)
{
  gmix[i].miu1=miu1[i];
  gmix[i].miu2=miu2[i];
  gmix[i].var1=var1[i];
  gmix[i].var2=var2[i];
  gmix[i].covar=covar[i];
  gmix[i].w=weight[i];
  gmix[i].ptr=new valtype[NofLocations];
}

std::cout<<"iteration and converge diff: ";

//std::vector<int>gmixIndexV(gmix.size());
//std::vector<int>locationIndexV(NofLocations);

// int icurrentVec[NofCore], *iCurrent=icurrentVec,
//   iBoundSize=NofCore+1,  iBoundsV[iBoundSize], *iStart=iBoundsV;
//
// signed char haultedV[NofCore], *haulted=haultedV;
//std::fill(haulted, haulted+NofCore, 0);

// volatile bool someoneQuit=0;

// std::vector<G2>gmix2=gmix;
// {
// for(int i=0,iend=gmix2.size();i!=iend;++i)
// {
//   gmix2[i].ptr=new valtype[NofLocations];
//   std::copy(gmix[i].ptr, gmix[i].ptr+NofLocations, gmix2[i].ptr);
// }
// }

dynamicTasking dT;

for(int k=0;k!=maxit;++k)
{
  //double converageDiff=0;


  // for(int i=0,iend=gmix2.size();i!=iend;++i)// evaluate densities with single core
  // {
  //   gmix2[i].evalV(&*Long.begin(), &*Lat.begin(), NofLocations, gmix2[i].ptr);
  // }




  int gmixSize=gmix.size();
  //
  double*LongBegin=&*Long.begin(), *LatBegin=&*Lat.begin();
  G2*gmixSt=&gmix.front();
  //
  // gmixIndexV.resize(gmix.size());
  // for(int i=0,iend=gmix.size();i!=iend;++i)gmixIndexV[i]=i;
  //
  // int avgtask=gmix.size()/NofCore;
  // for(int i=0,iend=NofCore;i!=iend;++i)
  // {
  //   iStart[i]=avgtask*i;
  //   iCurrent[i]=iStart[i];
  // }
  // iStart[NofCore]=gmix.size();
  //
  // std::fill(haulted, haulted+NofCore, 0);

  // firstForEvalDensities T(gmixSt, gmixSize, LongBegin, LatBegin, NofLocations,
  //                         NofCore, someoneQuit, gmixIndexV, iStart, iCurrent, haulted);

  //dT.reset(NofCore, gmix.size());
  //std::cout<<1.1<<std::endl;
  dT.reset(gmix.size());
  firstForEvalDensities T(gmixSt, LongBegin, LatBegin, NofLocations, ellipseAxisRatioThreshold, dT);

  parallelFor(0, NofCore, T);
  //std::cout<<1.2<<std::endl;

  // compare multicore and single core result
  // {
  // valtype tmpSS=0;
  // for(int i=0,iend=gmix.size();i!=iend;++i)
  // {
  //   valtype tmpS=0;
  //   for(int j=0,jend=NofLocations;j!=jend;++j)
  //   {
  //     valtype tmp=std::abs(gmix[i].ptr[j]-gmix2[i].ptr[j]);
  //     if(tmp!=0)
  //     {
  //       std::cout<<i<<" "<<j<<" "<<"diff="<<tmp<<" "<<"singleCore="<<gmix2[i].ptr[j]
  //                        <<" "<<"multicore="<<gmix[i].ptr[j]<<std::endl;
  //       std::cout<<"parameterDiff="<<std::abs(gmix[i].w-gmix2[i].w)+std::abs(gmix[i].miu1-gmix2[i].miu1)+
  //         std::abs(gmix[i].miu2-gmix2[i].miu2)+std::abs(gmix[i].covar-gmix2[i].covar)+
  //         std::abs(gmix[i].var1-gmix2[i].var1)+std::abs(gmix[i].var2-gmix2[i].var2)<<std::endl<<std::endl;
  //     }
  //     tmpS+=tmp;
  //   }
  //   tmpSS+=tmpS;
  // }
  // std::cout<<"initial density difference="<<tmpSS<<std::endl;
  // }





  // for(int j=0;j!=NofLocations;++j)
  // {
  //   valtype tmpS=0;
  //   for(int i=0,iend=gmix2.size();i!=iend;++i)
  //   {
  //     tmpS+=gmix2[i].ptr[j];
  //   }
  //   //std::cout<<"tmpS=="<<tmpS<<std::endl;
  //   for(int i=0,iend=gmix2.size();i!=iend;++i)
  //   {
  //     gmix2[i].ptr[j]/=tmpS;
  //   }
  // }

  //dT.reset(NofCore, NofLocations);
  dT.reset(NofLocations);

  //dynamicTasking dT2(NofCore, NofLocations);
  //struct secondForEvalWeightM T2(gmixSt, gmixSize, NofLocations, NofCore, dT);
  struct secondForEvalWeightM T2(gmixSt, gmixSize, dT);
  parallelFor(0, NofCore, T2);


  // {
  // valtype tmpSS=0;
  // for(int i=0,iend=gmix.size();i!=iend;++i)
  // {
  //   valtype tmpS=0;
  //   for(int j=0,jend=NofLocations;j!=jend;++j)
  //   {
  //     valtype tmp=std::abs(gmix[i].ptr[j]-gmix2[i].ptr[j]);
  //     //if(tmp!=0)std::cout<<i<<" "<<j<<std::endl;
  //     tmpS+=tmp;
  //   }
  //   tmpSS+=tmpS;
  // }
  // std::cout<<"normalized density difference="<<tmpSS<<std::endl;
  // }




  // for(int i=gmix2.size()-1;i>=0;--i)
  // {
  //   valtype miu1r=gmix2[i].miu1, miu2r=gmix2[i].miu2, var1r=gmix2[i].var1,
  //     var2r=gmix2[i].var2, covarr=gmix2[i].covar, wr=gmix2[i].w;
  //
  //   double*lossi=&*val.begin();
  //   for(int j=0;j!=NofLocations;++j,++lossi)
  //     gmix2[i].ptr[j]*=*lossi;
  //
  //   valtype Nk=std::accumulate(gmix2[i].ptr, gmix2[i].ptr+NofLocations, 0.0);
  //   gmix2[i].w=Nk/sumLoss;
  //
  //   if(gmix2[i].w<weightEPS)
  //   {
  //     //gmix2.erase(gmix2.begin()+i);
  //     continue;
  //   }
  //
  //   gmix2[i].miu1=0;
  //   gmix2[i].miu2=0;
  //   gmix2[i].var1=0;
  //   gmix2[i].var2=0;
  //   gmix2[i].covar=0;
  //   for(int j=0;j!=NofLocations;++j)
  //   {
  //     gmix2[i].miu1+=Long[j]*gmix2[i].ptr[j];
  //     gmix2[i].miu2+=Lat[j]*gmix2[i].ptr[j];
  //   }
  //   gmix2[i].miu1/=Nk;
  //   gmix2[i].miu2/=Nk;
  //
  //   for(int j=0;j!=NofLocations;++j)
  //   {
  //     gmix2[i].var1+=(Long[j]-gmix2[i].miu1)*(Long[j]-gmix2[i].miu1)*gmix2[i].ptr[j];
  //     gmix2[i].var2+=(Lat[j]-gmix2[i].miu2)*(Lat[j]-gmix2[i].miu2)*gmix2[i].ptr[j];
  //     gmix2[i].covar+=(Long[j]-gmix2[i].miu1)*(Lat[j]-gmix2[i].miu2)*gmix2[i].ptr[j];
  //   }
  //
  //   gmix2[i].var1/=Nk;
  //   gmix2[i].var2/=Nk;
  //   gmix2[i].covar/=Nk;
  // //
  // //   converageDiff+=std::abs(gmix2[i].var1-var1r)+std::abs(gmix2[i].var2-var2r)
  // //     +std::abs(gmix2[i].covar-covarr)+std::abs(gmix2[i].miu1-miu1r)+
  // //       std::abs(gmix2[i].miu2-miu2r)+std::abs(gmix2[i].w-wr);
  //  }


  double*valBegin=&*val.begin();
  valtype diffV[NofCore], *diff=diffV;
  std::fill(diff, diff+NofCore, 0);

  //dT.reset(NofCore, gmix.size());
  dT.reset(gmix.size());
  thirdForEvalGaussians T3(gmix, NofLocations, //NofCore,
                           valBegin, sumLoss,
                           weightEPS, LongBegin, LatBegin, diff,
                           convergeRelativeDiff, dT);

  parallelFor(0, NofCore, T3);


  for(int u=gmix.size()-1;u>=0;--u)
  {
    if(gmix[u].w<weightEPS)
    {
      delete[] gmix[u].ptr;
      gmix.erase(gmix.begin()+u);
    }
  }

  // {
  // valtype tmpS=0;
  // for(int i=0,iend=gmix.size();i!=iend;++i)
  // {
  //   //valtype tmpS=0;
  //   tmpS+=std::abs(gmix[i].w-gmix2[i].w)+std::abs(gmix[i].miu1-gmix2[i].miu1)+
  //     std::abs(gmix[i].miu2-gmix2[i].miu2)+
  //     std::abs(gmix[i].var1-gmix2[i].var1)+std::abs(gmix[i].covar-gmix2[i].covar)+
  //     std::abs(gmix[i].var2-gmix2[i].var2);
  // }
  // std::cout<<"parameter difference sum="<<tmpS<<std::endl;
  // }


  valtype avgDiff=std::accumulate(diffV, diffV+NofCore, 0.0)/(gmix.size()*6);
  std::cout<<k<<", "<<avgDiff<<"  ";

  if(avgDiff<convergeEPS)
  {
    //std::cout<<"meant to end\n";
    bool narrowComponentExist=0;
    for(int u=0,uend=gmix.size();u!=uend;++u)
    {
      if(gmix[u].eigenRatioApp()>ellipseAxisRatioThreshold)
      {
        std::cout<<gmix[u].eigenRatioApp()<<std::endl;
        narrowComponentExist=1;
        break;
      }
    }
    if(!narrowComponentExist)break;
  }

//
//
//
//   valtype avgDiff=std::accumulate(diffV, diffV+NofCore, 0.0)/(gmix.size()*6);
//   if(k%showProgress==0)
//   {
//     std::cout<<k<<", "<<avgDiff<<"   ";
//   }
//   if(avgDiff<convergeEPS)break;
}


List rst(gmix.size());
for(int i=0,iend=rst.size();i!=iend;++i)
{
  NumericVector result(6);
  result[0]=gmix[i].w;
  result[1]=gmix[i].miu1;
  result[2]=gmix[i].miu2;
  result[3]=gmix[i].var1;
  result[4]=gmix[i].covar;
  result[5]=gmix[i].var2;
  rst[i]=result;
  delete[] gmix[i].ptr;
}


// {
// for(int i=0,iend=rst.size();i!=iend;++i)
// {
//   delete[] gmix2[i].ptr;
// }
// }



return rst;
}












// [[Rcpp::export]]
List gm2d(NumericVector Long, NumericVector Lat, NumericVector val, NumericVector weight,
              NumericVector miu1, NumericVector miu2, NumericVector var1, NumericVector var2,
              NumericVector covar, double weightEPS, double convergeEPS,
              int maxit, int showProgress, bool useRelativeDiff, double ellipseAxisRatioThreshold=14.2){

int intialG=miu1.size();
//std::cout<<"intialG="<<intialG<<std::endl;
valtype sumLoss=std::accumulate(val.begin(),val.end(),0.0);
int NofLocations=Long.size();

std::vector<G2>gmix(intialG);
//std::vector<std::vector<valtype> >wM(intialG, std::vector<valtype>(NofLocations));


for(int i=0,iend=gmix.size();i!=iend;++i)
{
  gmix[i].miu1=miu1[i];
  gmix[i].miu2=miu2[i];
  gmix[i].var1=var1[i];
  gmix[i].var2=var2[i];
  gmix[i].covar=covar[i];
  gmix[i].w=weight[i];
  //gmix[i].ptr=&wM[i].front();
  gmix[i].ptr=new valtype[NofLocations];
}

std::cout<<"iteration and converge diff: ";
for(int k=0;k!=maxit;++k)
{
  double converageDiff=0;

  for(int i=0,iend=gmix.size();i!=iend;++i)
  {
    //std::cout<<"i=="<<i<<std::endl;
    gmix[i].evalV(&*Long.begin(), &*Lat.begin(), NofLocations, gmix[i].ptr, ellipseAxisRatioThreshold);
  }

  for(int j=0;j!=NofLocations;++j)
  {
    valtype tmpS=0;
    for(int i=0,iend=gmix.size();i!=iend;++i)
    {
      tmpS+=gmix[i].ptr[j];
    }
    //std::cout<<"tmpS=="<<tmpS<<std::endl;
    for(int i=0,iend=gmix.size();i!=iend;++i)
    {
      gmix[i].ptr[j]/=tmpS;
      if(!std::isfinite(gmix[i].ptr[j]))gmix[i].ptr[j]=0;
    }
  }

  // output the membership weight matrix
  // for(int j=0;j!=NofLocations;++j)
  // {
  //   for(int i=0,iend=gmix.size();i!=iend;++i)
  //     std::cout<<gmix[i].ptr[j]<<" ";
  //   std::cout<<std::endl;
  // }


  int progress=0;
  for(int i=gmix.size()-1;i>=0;--i)
  {

    G2 tmpReserve=gmix[i];

    // valtype miu1r=gmix[i].miu1, miu2r=gmix[i].miu2, var1r=gmix[i].var1,
    //   var2r=gmix[i].var2, covarr=gmix[i].covar, wr=gmix[i].w;

    double*lossi=&*val.begin();
    for(int j=0;j!=NofLocations;++j,++lossi)
      gmix[i].ptr[j]*=*lossi;

    valtype Nk=std::accumulate(gmix[i].ptr, gmix[i].ptr+NofLocations, 0.0);
    gmix[i].w=Nk/sumLoss;

    if(gmix[i].w<weightEPS)
    {
      delete[] gmix[i].ptr;
      gmix.erase(gmix.begin()+i);
      continue;
    }

    gmix[i].miu1=0;
    gmix[i].miu2=0;
    gmix[i].var1=0;
    gmix[i].var2=0;
    gmix[i].covar=0;
    for(int j=0;j!=NofLocations;++j)
    {
      gmix[i].miu1+=Long[j]*gmix[i].ptr[j];
      gmix[i].miu2+=Lat[j]*gmix[i].ptr[j];
    }
    gmix[i].miu1/=Nk;
    gmix[i].miu2/=Nk;

    for(int j=0;j!=NofLocations;++j)
    {
      gmix[i].var1+=(Long[j]-gmix[i].miu1)*(Long[j]-gmix[i].miu1)*gmix[i].ptr[j];
      gmix[i].var2+=(Lat[j]-gmix[i].miu2)*(Lat[j]-gmix[i].miu2)*gmix[i].ptr[j];
      gmix[i].covar+=(Long[j]-gmix[i].miu1)*(Lat[j]-gmix[i].miu2)*gmix[i].ptr[j];
    }

    gmix[i].var1/=Nk;
    gmix[i].var2/=Nk;
    gmix[i].covar/=Nk;

    // converageDiff+=std::abs(gmix[i].var1-var1r)+std::abs(gmix[i].var2-var2r)
    //   +std::abs(gmix[i].covar-covarr)+std::abs(gmix[i].miu1-miu1r)+
    //     std::abs(gmix[i].miu2-miu2r)+std::abs(gmix[i].w-wr);

    converageDiff+=evalDiff(tmpReserve, gmix[i], useRelativeDiff);
    ++progress;
  }

  if(k%showProgress==0)
  {
    std::cout<<k<<", "<<converageDiff/(progress*6)<<"   ";
  }

  if(converageDiff<convergeEPS*(progress*6))
  {
    //std::cout<<"meant to end\n";
    bool narrowComponentExist=0;
    for(int u=0,uend=gmix.size();u!=uend;++u)
    {
      if(gmix[u].eigenRatioApp()>ellipseAxisRatioThreshold)
      {
        std::cout<<gmix[u].eigenRatioApp()<<std::endl;
        narrowComponentExist=1;
        break;
      }
    }
    if(!narrowComponentExist)break;
  }
}


List rst(gmix.size());
for(int i=0,iend=rst.size();i!=iend;++i)
{
  NumericVector result(6);
  result[0]=gmix[i].w;
  result[1]=gmix[i].miu1;
  result[2]=gmix[i].miu2;
  result[3]=gmix[i].var1;
  result[4]=gmix[i].covar;
  result[5]=gmix[i].var2;
  rst[i]=result;
  delete[] gmix[i].ptr;
}

return rst;
}











// [[Rcpp::export]]
List gm2dComponentWise(NumericVector Long, NumericVector Lat, NumericVector val, NumericVector weight,
              NumericVector miu1, NumericVector miu2, NumericVector var1, NumericVector var2,
              NumericVector covar, double weightEPS, double convergeEPS,
              int maxit, int showProgress, bool useRelativeDiff, double ellipseAxisRatioThreshold=14.2){

int intialG=miu1.size();
//std::cout<<"intialG="<<intialG<<std::endl;
valtype sumLoss=std::accumulate(val.begin(),val.end(),0.0);
int NofLocations=Long.size();

std::vector<G2>gmix(intialG);
//std::vector<std::vector<valtype> >wM(intialG, std::vector<valtype>(NofLocations));


for(int i=0,iend=gmix.size();i!=iend;++i)
{
  gmix[i].miu1=miu1[i];
  gmix[i].miu2=miu2[i];
  gmix[i].var1=var1[i];
  gmix[i].var2=var2[i];
  gmix[i].covar=covar[i];
  gmix[i].w=weight[i];
  //gmix[i].ptr=&wM[i].front();
  gmix[i].ptr=new valtype[NofLocations];
}


for(int i=0,iend=gmix.size();i!=iend;++i)
{
    //std::cout<<"i=="<<i<<std::endl;
  gmix[i].evalV(&*Long.begin(), &*Lat.begin(), NofLocations, gmix[i].ptr,
                  ellipseAxisRatioThreshold);
}






std::cout<<"iteration and converge diff: ";
for(int k=0;k!=maxit;++k)
{
  double converageDiff=0;

  // for(int i=0,iend=gmix.size();i!=iend;++i)
  // {
  //   //std::cout<<"i=="<<i<<std::endl;
  //   gmix[i].evalV(&*Long.begin(), &*Lat.begin(), NofLocations, gmix[i].ptr,
  //                 ellipseAxisRatioThreshold);
  // }
  //
  // for(int j=0;j!=NofLocations;++j)
  // {
  //   valtype tmpS=0;
  //   for(int i=0,iend=gmix.size();i!=iend;++i)
  //   {
  //     tmpS+=gmix[i].ptr[j];
  //   }
  //   //std::cout<<"tmpS=="<<tmpS<<std::endl;
  //   for(int i=0,iend=gmix.size();i!=iend;++i)
  //   {
  //     gmix[i].ptr[j]/=tmpS;
  //     if(!std::isfinite(gmix[i].ptr[j]))gmix[i].ptr[j]=0;
  //   }
  // }

  // output the membership weight matrix
  // for(int j=0;j!=NofLocations;++j)
  // {
  //   for(int i=0,iend=gmix.size();i!=iend;++i)
  //     std::cout<<gmix[i].ptr[j]<<" ";
  //   std::cout<<std::endl;
  // }


  int progress=0;
  //bool nextIter=0;

  for(int i=gmix.size()-1;i>=0;--i)
  {
    // valtype miu1r=gmix[i].miu1, miu2r=gmix[i].miu2, var1r=gmix[i].var1,
    //   var2r=gmix[i].var2, covarr=gmix[i].covar, wr=gmix[i].w;

    for(int j=0;j!=NofLocations;++j)// calculate the membership matrix
    {
      valtype tmpS=0;
      for(int u=0,uend=gmix.size();u!=uend;++u)
        tmpS+=gmix[u].ptr[j];
      for(int u=0,uend=gmix.size();u!=uend;++u)
      {
        gmix[u].ptr[j]/=tmpS;
        if(!std::isfinite(gmix[u].ptr[j]))gmix[u].ptr[j]=0;
      }
    }


    G2 tmpReserve=gmix[i];

    double*lossi=&*val.begin();
    for(int j=0;j!=NofLocations;++j,++lossi)
      gmix[i].ptr[j]*=*lossi;

    valtype Nk=std::accumulate(gmix[i].ptr, gmix[i].ptr+NofLocations, 0.0);
    gmix[i].w=Nk/sumLoss;

    if(gmix[i].w<weightEPS)
    {
      delete[] gmix[i].ptr;
      gmix.erase(gmix.begin()+i);
      continue;
    }

    gmix[i].miu1=0;
    gmix[i].miu2=0;
    gmix[i].var1=0;
    gmix[i].var2=0;
    gmix[i].covar=0;
    for(int j=0;j!=NofLocations;++j)
    {
      gmix[i].miu1+=Long[j]*gmix[i].ptr[j];
      gmix[i].miu2+=Lat[j]*gmix[i].ptr[j];
    }
    gmix[i].miu1/=Nk;
    gmix[i].miu2/=Nk;

    for(int j=0;j!=NofLocations;++j)
    {
      gmix[i].var1+=(Long[j]-gmix[i].miu1)*(Long[j]-gmix[i].miu1)*gmix[i].ptr[j];
      gmix[i].var2+=(Lat[j]-gmix[i].miu2)*(Lat[j]-gmix[i].miu2)*gmix[i].ptr[j];
      gmix[i].covar+=(Long[j]-gmix[i].miu1)*(Lat[j]-gmix[i].miu2)*gmix[i].ptr[j];
    }

    gmix[i].var1/=Nk;
    gmix[i].var2/=Nk;
    gmix[i].covar/=Nk;

    for(int u=0,uend=gmix.size();u!=uend;++u)// recompute the densities
    {
      gmix[u].evalV(&*Long.begin(), &*Lat.begin(), NofLocations, gmix[u].ptr,
          ellipseAxisRatioThreshold);
    }

    converageDiff+=evalDiff(tmpReserve, gmix[i], useRelativeDiff);
    ++progress;
  }


  if(k%showProgress==0)std::cout<<k<<", "<<converageDiff/(progress*6)<<"   ";
  if(converageDiff<convergeEPS*(progress*6))break;
}


List rst(gmix.size());
for(int i=0,iend=rst.size();i!=iend;++i)
{
  NumericVector result(6);
  result[0]=gmix[i].w;
  result[1]=gmix[i].miu1;
  result[2]=gmix[i].miu2;
  result[3]=gmix[i].var1;
  result[4]=gmix[i].covar;
  result[5]=gmix[i].var2;
  rst[i]=result;
  delete[] gmix[i].ptr;
}

return rst;
}


























































valtype kernelcovar(G2&X, G2&Y){
valtype sigmaSum[4];
sigmaSum[0]=X.var1+Y.var1;
sigmaSum[1]=X.covar+Y.covar;
sigmaSum[2]=sigmaSum[1];
sigmaSum[3]=X.var2+Y.var2;
valtype m=sigmaSum[0]*sigmaSum[3]-sigmaSum[1]*sigmaSum[2];
valtype denominator=std::sqrt(m)*M_PI*2;
std::swap(sigmaSum[0],sigmaSum[3]);
sigmaSum[0]/=m;
sigmaSum[1]/=-m;
sigmaSum[2]/=-m;
sigmaSum[3]/=m;
valtype rst=std::exp((sigmaSum[0]*(X.miu1-Y.miu1)*(X.miu1-Y.miu1)+(sigmaSum[1]+sigmaSum[2])*
  (X.miu1-Y.miu1)*(X.miu2-Y.miu2)+sigmaSum[3]*(X.miu2-Y.miu2)*(X.miu2-Y.miu2))/(-2))
  /denominator;
if(!std::isfinite(rst))rst=0;
return rst;
}










void kernelcovar(std::vector<G2>&x, std::vector<valtype>&kernelConvCovar){
unsigned xsize=x.size();
kernelConvCovar.resize(xsize*xsize);
for(unsigned i=0,iend=xsize;i!=iend;++i)
{
  for(unsigned j=0,jend=xsize;j!=jend;++j)
    kernelConvCovar[i*xsize+j]=kernelcovar(x[i],x[j]);
}
}








void locationKernelIndWeight(double*x, double*y, unsigned NofLocation,
                 std::vector<std::vector<unsigned> >&ind,
                 std::vector<std::vector<valtype> >&w,
                 std::vector<G2>&kernel,
                 double weightEPS, bool useMahalanobisD){

ind.resize(NofLocation);
w.resize(ind.size());

for(unsigned i=0,iend=ind.size();i!=iend;++i)
{
  //if(i%10000==0)std::cout<<i<<" ";
  std::vector<valtype>tmpw(kernel.size());
  std::vector<unsigned>tmpind(kernel.size());
  //ind[i].resize(kernel.size());
  //valtype tmpwS=0;
  for(unsigned j=0,jend=kernel.size();j!=jend;++j)
  {
    tmpind[j]=j;
    if(!useMahalanobisD)
      tmpw[j]=kernel[j].evalNoWeight(x[i],y[i]);
    else tmpw[j]=kernel[j].evalMahalanobisD(x[i],y[i]);
    //tmpwS+=tmpw[j];
  }
  valtype tmpwS=0;

  // here comes the variation!
  // (1):
  for(unsigned j=0,jend=kernel.size();j!=jend;++j)tmpwS+=tmpw[j];
  for(unsigned j=0,jend=kernel.size();j!=jend;++j)tmpw[j]/=tmpwS;

  // // (2):
  // for(unsigned j=0,jend=kernel.size();j!=jend;++j)tmpwS+=tmpw[j]*std::abs(kernel[j].w);
  // for(unsigned j=0,jend=kernel.size();j!=jend;++j)tmpw[j]*=std::abs(kernel[j].w)/tmpwS;

  // // (3):
  // for(unsigned j=0,jend=kernel.size();j!=jend;++j)tmpwS+=tmpw[j]*std::abs(kernel[j].w);
  // for(unsigned j=0,jend=kernel.size();j!=jend;++j)tmpw[j]*=kernel[j].w*kernel[j].w/tmpwS;

  tmpwS=0;
  for(int j=kernel.size()-1;j>=0;--j)
  {
    if(tmpw[j]<weightEPS)
    {
      tmpw.erase(tmpw.begin()+j);
      tmpind.erase(tmpind.begin()+j);
    }
    else tmpwS+=tmpw[j];
  }
  for(unsigned j=0,jend=tmpw.size();j!=jend;++j)tmpw[j]/=tmpwS;
  ind[i].assign(tmpind.begin(), tmpind.end());
  w[i].assign(tmpw.begin(),tmpw.end());
}
}










void computelocationKernelSD(double*x, double*y, unsigned NofLocation,
                 std::vector<G2>&kernel, double weightEPS,
                 std::vector<valtype>&locationKernelSD,
                 std::vector<valtype>&kernelConvCovar,
                 std::vector<std::vector<unsigned> >&IND,
                 std::vector<std::vector<valtype> >&W){

locationKernelSD.resize(NofLocation);

//std::vector<std::vector<unsigned> >IND;
//std::vector<std::vector<valtype> >W;
//locationKernelIndWeight(x, y, NofLocation, IND, W, kernel, weightEPS);

unsigned kernelSize=kernel.size();

for(unsigned i=0,iend=NofLocation;i!=iend;++i)
{
  unsigned*ind=&IND[i][0];
  valtype*w=&W[i][0];
  unsigned siz=IND[i].size();
  valtype covar=0;
  for(unsigned p=0;p!=siz;++p)
  {
    unsigned tmpInd=ind[p]*kernelSize;
    valtype covar2=0;
    for(unsigned q=0;q!=siz;++q)covar2+=kernelConvCovar[tmpInd+ind[q]]*w[p]*w[q];
    covar+=covar2;
  }
  locationKernelSD[i]=std::sqrt(covar);
}
}











valtype totalVar(double*x, double*y, unsigned NofLocation,
                 std::vector<G2>&kernel, double weightEPS,
                 double*sd, valtype sumVar, bool outputCorM, double*corM,
                 bool useMahalanobisD){

std::vector<std::vector<unsigned> >ind;
std::vector<std::vector<valtype> >w;
locationKernelIndWeight(x, y, NofLocation, ind, w, kernel, weightEPS, useMahalanobisD);

unsigned kernelSize=kernel.size();
std::vector<valtype>kernelConvCovar;

kernelcovar(kernel, kernelConvCovar);

std::vector<valtype>locationKernelSD;

computelocationKernelSD(x, y, NofLocation, kernel, weightEPS,
                        locationKernelSD, kernelConvCovar, ind, w);

valtype covar=0;
for(unsigned i=0,iend=NofLocation-1;i!=iend;++i)
{
  //if(i%256==0)std::cout<<i<<" ";
  valtype covar2=0;
  for(unsigned j=i+1,jend=NofLocation;j!=jend;++j)
  {
    unsigned*ind1=&ind[i][0], *ind2=&ind[j][0];
    valtype*w1=&w[i][0], *w2=&w[j][0];
    unsigned siz1=ind[i].size(), siz2=ind[j].size();
    valtype covar3=0;
    for(unsigned p=0;p!=siz1;++p)
    {
      unsigned tmpInd=ind1[p]*kernelSize;
      valtype covar4=0;
      for(unsigned q=0;q!=siz2;++q)
      {
        covar4+=kernelConvCovar[tmpInd+ind2[q]]*w1[p]*w2[q];
      }
      covar3+=covar4;
    }
    //covar3 is the covariance of i and j by Gaussian Process
    valtype corr=covar3/(locationKernelSD[i]*locationKernelSD[j]);
    covar2+=corr*(sd[i]*sd[j]);
    if(outputCorM)
    {
      corM[NofLocation*i+j]=corr;
      corM[NofLocation*j+i]=corr;
    }
  }
  covar+=covar2;
}
if(outputCorM)
{
  for(unsigned i=0;i!=NofLocation;++i)corM[i*NofLocation+i]=1;
}
return sumVar+covar+covar;
}










// [[Rcpp::export]]
double GMKtotalVarCpp(List location, List GMkernel, NumericVector locationVar,
                      NumericMatrix corM, double weightEPS=1e-12, bool outputCorM=0,
                      bool useMahalanobisD=0){

std::vector<G2>kernel(GMkernel.size());
for(int i=0,iend=kernel.size();i!=iend;++i)
{
  NumericVector tmp=GMkernel[i];
  kernel[i].w=tmp[0];
  kernel[i].miu1=tmp[1];
  kernel[i].miu2=tmp[2];
  kernel[i].var1=tmp[3];
  kernel[i].covar=tmp[4];
  kernel[i].var2=tmp[5];
}
NumericVector x=location[0], y=location[1];

std::vector<valtype>locationSD(locationVar.size());
for(unsigned i=0,iend=locationSD.size();i!=iend;++i)
  locationSD[i]=std::sqrt(locationVar[i]);

return totalVar(&x[0], &y[0], locationVar.size(), kernel, weightEPS,
                &locationSD.front(), std::accumulate(locationVar.begin(),
                locationVar.end(),0.0), outputCorM, &*corM.begin(), useMahalanobisD);

// return totalVarAndOutputCorM(&x[0], &y[0], locationVar.size(),
//    kernel, weightEPS, &locationSD.front(),
//    std::accumulate(locationVar.begin(),locationVar.end(),0.0), &*corM.begin());
}









struct NS: public Worker{
double*&x;
double*&y;
unsigned& NofLocation;
std::vector<G2>&kernel;
double& weightEPS;
double*&sd;
valtype& sumVar;
bool& outputCorM;
double*&corM;
unsigned*&index;
std::vector<std::vector<unsigned> >&ind;
std::vector<std::vector<valtype> >&w;
std::vector<valtype>&kernelConvCovar;
std::vector<valtype>&locationKernelSD;
valtype*&covarSum;
unsigned*&NofColumnProcessed;
unsigned&NofCore;

//Here is how to initialize class member of references
NS(double*&x,
double*&y,
unsigned& NofLocation,
std::vector<G2>&kernel,
double& weightEPS,
double*&sd,
valtype& sumVar,
bool& outputCorM,
double*&corM,
unsigned*&index,
std::vector<std::vector<unsigned> >&ind,
std::vector<std::vector<valtype> >&w,
std::vector<valtype>&kernelConvCovar,
std::vector<valtype>&locationKernelSD,
valtype*&covarSum,
unsigned*&NofColumnProcessed,
unsigned&NofCore):

x(x),
y(y),
NofLocation(NofLocation),
kernel(kernel),
weightEPS(weightEPS),
sd(sd),
sumVar(sumVar),
outputCorM(outputCorM),
corM(corM),
index(index),
ind(ind),
w(w),
kernelConvCovar(kernelConvCovar),
locationKernelSD(locationKernelSD),
covarSum(covarSum),
NofColumnProcessed(NofColumnProcessed),
NofCore(NofCore){}

//when writing operator() function, it must use std::size_t!! Anything else won't pass!!
void operator()(std::size_t st, std::size_t end){
for(std::size_t I=st;I!=end;++I)
{
unsigned i, iend=index[I];
if(I==0)i=0;
else i=index[I-1];

unsigned kernelSize=kernel.size();
valtype covar=0;

std::clock_t timestart=std::clock();

for(;i!=iend;++i,++NofColumnProcessed[I])
{
  std::clock_t timeend=std::clock();
  if((timeend-timestart)/(double)CLOCKS_PER_SEC>300)
  {
    std::cout<<-std::accumulate(NofColumnProcessed, NofColumnProcessed+NofCore, 0);
    timestart=timeend;
  }

  valtype covar2=0;
  for(unsigned j=i+1,jend=NofLocation;j!=jend;++j)
  {
    unsigned*ind1=&ind[i][0], *ind2=&ind[j][0];
    valtype*w1=&w[i][0], *w2=&w[j][0];
    unsigned siz1=ind[i].size(), siz2=ind[j].size();
    valtype covar3=0;
    for(unsigned p=0;p!=siz1;++p)
    {
      unsigned tmpInd=ind1[p]*kernelSize;
      valtype covar4=0;
      for(unsigned q=0;q!=siz2;++q)
      {
        covar4+=kernelConvCovar[tmpInd+ind2[q]]*w1[p]*w2[q];
      }
      covar3+=covar4;
    }
    //covar3 is the covariance of i and j by Gaussian Process
    valtype corr=covar3/(locationKernelSD[i]*locationKernelSD[j]);
    covar2+=corr*(sd[i]*sd[j]);
    if(outputCorM)
    {
      corM[NofLocation*i+j]=corr;
      corM[NofLocation*j+i]=corr;
    }
  }
  covar+=covar2;
}
covarSum[I]=covar;
}
}
};






struct NSfightForColumn: public Worker{
double*&x;
double*&y;
unsigned& NofLocation;
std::vector<G2>&kernel;
double& weightEPS;
double*&sd;
valtype& sumVar;
bool& outputCorM;
double*&corM;
//unsigned*&index;
std::vector<std::vector<unsigned> >&ind;
std::vector<std::vector<valtype> >&w;
std::vector<valtype>&kernelConvCovar;
std::vector<valtype>&locationKernelSD;
valtype*&covarSum;
unsigned*&NofColumnProcessed;
unsigned&NofCore;
//unsigned char*&ifcolprocessed;
volatile unsigned&availableCol;

//Here is how to initialize class member of references
NSfightForColumn(double*&x,
double*&y,
unsigned& NofLocation,
std::vector<G2>&kernel,
double& weightEPS,
double*&sd,
valtype& sumVar,
bool& outputCorM,
double*&corM,
//unsigned*&index,
std::vector<std::vector<unsigned> >&ind,
std::vector<std::vector<valtype> >&w,
std::vector<valtype>&kernelConvCovar,
std::vector<valtype>&locationKernelSD,
valtype*&covarSum,
unsigned*&NofColumnProcessed,
unsigned&NofCore,
//unsigned char*&ifcolprocessed
volatile unsigned&availableCol):

x(x),
y(y),
NofLocation(NofLocation),
kernel(kernel),
weightEPS(weightEPS),
sd(sd),
sumVar(sumVar),
outputCorM(outputCorM),
corM(corM),
//index(index),
ind(ind),
w(w),
kernelConvCovar(kernelConvCovar),
locationKernelSD(locationKernelSD),
covarSum(covarSum),
NofColumnProcessed(NofColumnProcessed),
NofCore(NofCore),
//ifcolprocessed(ifcolprocessed)
availableCol(availableCol){}

//when writing operator() function, it must use std::size_t!! Anything else won't pass!!
void operator()(std::size_t st, std::size_t end){
for(std::size_t I=st;I!=end;++I)
{
  valtype covar=0;
  unsigned kernelSize=kernel.size();

for(;availableCol<NofLocation;)
{
  unsigned i=availableCol;
  ++availableCol;

  valtype covar2=0;
  for(unsigned j=i+1,jend=NofLocation;j!=jend;++j)
  {
    unsigned*ind1=&ind[i][0], *ind2=&ind[j][0];
    valtype*w1=&w[i][0], *w2=&w[j][0];
    unsigned siz1=ind[i].size(), siz2=ind[j].size();
    valtype covar3=0;
    for(unsigned p=0;p!=siz1;++p)
    {
      unsigned tmpInd=ind1[p]*kernelSize;
      valtype covar4=0;
      for(unsigned q=0;q!=siz2;++q)
      {
        covar4+=kernelConvCovar[tmpInd+ind2[q]]*w1[p]*w2[q];
      }
      covar3+=covar4;
    }
    //covar3 is the covariance of i and j by Gaussian Process
    valtype corr=covar3/(locationKernelSD[i]*locationKernelSD[j]);
    covar2+=corr*(sd[i]*sd[j]);
    if(outputCorM)
    {
      corM[NofLocation*i+j]=corr;
      corM[NofLocation*j+i]=corr;
    }
  }
  covar+=covar2;
}

covarSum[I]=covar;
}
}
};







struct NSredistribute: public Worker{
double*&x;
double*&y;
unsigned& NofLocation;
std::vector<G2>&kernel;
double& weightEPS;
double*&sd;
valtype& sumVar;
bool& outputCorM;
double*&corM;
unsigned*&indexV;
std::vector<std::vector<unsigned> >&ind;
std::vector<std::vector<valtype> >&w;
std::vector<valtype>&kernelConvCovar;
std::vector<valtype>&locationKernelSD;
valtype*&covarSum;
unsigned*&NofColumnProcessed;
unsigned&NofCore;
volatile bool&someoneQuit;
unsigned&NcolLeft;
unsigned*&leftIbegin;

//Here is how to initialize class member of references
NSredistribute(double*&x,
double*&y,
unsigned& NofLocation,
std::vector<G2>&kernel,
double& weightEPS,
double*&sd,
valtype& sumVar,
bool& outputCorM,
double*&corM,
unsigned*&indexV,
std::vector<std::vector<unsigned> >&ind,
std::vector<std::vector<valtype> >&w,
std::vector<valtype>&kernelConvCovar,
std::vector<valtype>&locationKernelSD,
valtype*&covarSum,
unsigned*&NofColumnProcessed,
unsigned&NofCore,
volatile bool&someoneQuit,
unsigned&NcolLeft,
unsigned*&leftIbegin):

x(x),
y(y),
NofLocation(NofLocation),
kernel(kernel),
weightEPS(weightEPS),
sd(sd),
sumVar(sumVar),
outputCorM(outputCorM),
corM(corM),
indexV(indexV),
ind(ind),
w(w),
kernelConvCovar(kernelConvCovar),
locationKernelSD(locationKernelSD),
covarSum(covarSum),
NofColumnProcessed(NofColumnProcessed),
NofCore(NofCore),
someoneQuit(someoneQuit),
NcolLeft(NcolLeft),
leftIbegin(leftIbegin){}

//when writing operator() function, it must use std::size_t!! Anything else won't pass!!
void operator()(std::size_t st, std::size_t end){
for(std::size_t I=st;I!=end;++I)
{
unsigned avgColForEachCore=NcolLeft/NofCore;
unsigned indexVi=I*avgColForEachCore, indexViend=indexVi+avgColForEachCore;

if(I==NofCore-1)indexViend=NcolLeft;

unsigned kernelSize=kernel.size();
valtype covar=0;

std::clock_t timestart=std::clock();

for(;;++indexVi,++NofColumnProcessed[I])
{
  if(someoneQuit||indexVi==indexViend)
  {
    leftIbegin[I]=indexVi;
    break;
  }

  unsigned i=indexV[indexVi];

  std::clock_t timeend=std::clock();
  if((timeend-timestart)/(double)CLOCKS_PER_SEC>300)
  {
    std::cout<<-std::accumulate(NofColumnProcessed, NofColumnProcessed+NofCore, 0);
    timestart=timeend;
  }

  valtype covar2=0;
  for(unsigned j=i+1,jend=NofLocation;j!=jend;++j)
  {
    unsigned*ind1=&ind[i][0], *ind2=&ind[j][0];
    valtype*w1=&w[i][0], *w2=&w[j][0];
    unsigned siz1=ind[i].size(), siz2=ind[j].size();
    valtype covar3=0;
    for(unsigned p=0;p!=siz1;++p)
    {
      unsigned tmpInd=ind1[p]*kernelSize;
      valtype covar4=0;
      for(unsigned q=0;q!=siz2;++q)
      {
        covar4+=kernelConvCovar[tmpInd+ind2[q]]*w1[p]*w2[q];
      }
      covar3+=covar4;
    }
    //covar3 is the covariance of i and j by Gaussian Process
    valtype corr=covar3/(locationKernelSD[i]*locationKernelSD[j]);
    covar2+=corr*(sd[i]*sd[j]);
    if(outputCorM)
    {
      corM[NofLocation*i+j]=corr;
      corM[NofLocation*j+i]=corr;
    }
  }
  covar+=covar2;

}

covarSum[I]+=covar;
someoneQuit=1;
}
}
};














valtype totalVarParallel(double*x, double*y, unsigned NofLocation,
                 std::vector<G2>&kernel, double weightEPS,
                 double*sd, valtype sumVar, bool outputCorM, double*corM,
                 unsigned NofCore, unsigned*index,// bool equalDivide,
                 std::vector<unsigned>&indexVec, bool useMahalanobisDforWeight){

std::vector<std::vector<unsigned> >ind;
std::vector<std::vector<valtype> >w;
locationKernelIndWeight(x, y, NofLocation, ind, w, kernel, weightEPS, useMahalanobisDforWeight);

//unsigned kernelSize=kernel.size();
std::vector<valtype>kernelConvCovar;

kernelcovar(kernel, kernelConvCovar);

std::vector<valtype>locationKernelSD;

computelocationKernelSD(x, y, NofLocation, kernel, weightEPS,
                        locationKernelSD, kernelConvCovar, ind, w);


std::vector<valtype>covarVec(NofCore,0);
valtype*covarSum=&covarVec.front();

std::vector<unsigned>NofColumnProcessedV(NofCore,0);
unsigned*NofColumnProcessed=&NofColumnProcessedV.front();

// if(equalDivide)
// {
//   NS T(x,y,NofLocation,kernel,weightEPS,sd,sumVar,outputCorM,corM,index,ind,w,
//     kernelConvCovar,locationKernelSD,covarSum,NofColumnProcessed,NofCore);
//   parallelFor(0, NofCore, T);
// }
// // else
// // {
// //   unsigned availableCol=0;
// //   NSfightForColumn T(x,y,NofLocation,kernel,weightEPS,sd,sumVar,outputCorM,corM,ind,w,
// //     kernelConvCovar,locationKernelSD,covarSum,NofColumnProcessed,NofCore,availableCol);
// //   parallelFor(0, NofCore, T);
// // }
// else
{

  unsigned c1[NofCore], *coreStart=c1;
  unsigned c2[NofCore], *leftIbegin=c2;

  while(true)
  {
    unsigned*indexV=&indexVec.front();
    volatile bool someoneQuit=0;

    unsigned NcolLeft=indexVec.size();
    //if(NcolLeft==0)break;

    if(NcolLeft<NofCore)NofCore=1;

    unsigned avgLoadEachCore=NcolLeft/NofCore;
    for(unsigned i=0;i!=NofCore;++i)coreStart[i]=i*avgLoadEachCore;

    //std::cout<<NofCore<<" "<<NcolLeft<<"\n";

    NSredistribute T(x,y,NofLocation,kernel,weightEPS,sd,sumVar,outputCorM,corM,indexV,ind,w,
      kernelConvCovar,locationKernelSD,covarSum,NofColumnProcessed,NofCore,someoneQuit,
        NcolLeft,leftIbegin);

    parallelFor(0, NofCore, T);

    //std::cout<<1.1<<std::endl;

    if(NofCore==1)break;

    for(int i=NofCore-1;i>=0;--i)
    {
      //std::cout<<"bounds=="<<coreStart[i]<<" "<<leftIbegin[i]<<std::endl;
      indexVec.erase(indexVec.begin()+coreStart[i], indexVec.begin()+leftIbegin[i]);
    }
    //std::cout<<1.2<<std::endl;

  }
}


valtype covar=accumulate(covarVec.begin(), covarVec.end(), 0.0);
if(outputCorM)
{
  for(unsigned i=0;i!=NofLocation;++i)corM[i*NofLocation+i]=1;
}
return sumVar+covar+covar;
}














valtype kernelcorr(G2&X, G2&Y){
//valtype coef=std::pow((X.var1*X.var2-X.covar*X.covar)*(Y.var1*Y.var2-Y.covar*Y.covar),0.25);
valtype sigmaSum[4];
sigmaSum[0]=(X.var1+Y.var1)/2;
sigmaSum[1]=(X.covar+Y.covar)/2;
sigmaSum[2]=sigmaSum[1];
sigmaSum[3]=(X.var2+Y.var2)/2;
valtype m=sigmaSum[0]*sigmaSum[3]-sigmaSum[1]*sigmaSum[2];
valtype coef=std::pow((X.var1*X.var2-X.covar*X.covar)*
  (Y.var1*Y.var2-Y.covar*Y.covar),1.0/4)/std::sqrt(m);
//valtype denominator=std::sqrt(m)*M_PI*2;
std::swap(sigmaSum[0],sigmaSum[3]);
sigmaSum[0]/=m;
sigmaSum[1]/=-m;
sigmaSum[2]/=-m;
sigmaSum[3]/=m;

valtype rst=coef*std::exp(-(sigmaSum[0]*(X.miu1-Y.miu1)*(X.miu1-Y.miu1)+
                     (sigmaSum[1]+sigmaSum[2])*(X.miu1-Y.miu1)*(X.miu2-Y.miu2)+
                     sigmaSum[3]*(X.miu2-Y.miu2)*(X.miu2-Y.miu2)));

if(!std::isfinite(rst))rst=0;
return rst;
}









struct NSredistributeLocovar: public Worker{
double*&x;
double*&y;
unsigned& NofLocation;
//std::vector<G2>&kernel;
//double& weightEPS;
double*&sd;
valtype& sumVar;
bool& outputCorM;
double*&corM;
unsigned*&indexV;
//std::vector<std::vector<unsigned> >&ind;
//std::vector<std::vector<valtype> >&w;
//std::vector<valtype>&kernelConvCovar;
//std::vector<valtype>&locationKernelSD;
std::vector<G2>&localGkernel;
valtype*&covarSum;
unsigned*&NofColumnProcessed;
unsigned&NofCore;
volatile bool&someoneQuit;
unsigned&NcolLeft;
unsigned*&leftIbegin;

//Here is how to initialize class member of references
NSredistributeLocovar(double*&x,
double*&y,
unsigned& NofLocation,
//std::vector<G2>&kernel,
//double& weightEPS,
double*&sd,
valtype& sumVar,
bool& outputCorM,
double*&corM,
unsigned*&indexV,
//std::vector<std::vector<unsigned> >&ind,
//std::vector<std::vector<valtype> >&w,
//std::vector<valtype>&kernelConvCovar,
//std::vector<valtype>&locationKernelSD,
std::vector<G2>&localGkernel,
valtype*&covarSum,
unsigned*&NofColumnProcessed,
unsigned&NofCore,
volatile bool&someoneQuit,
unsigned&NcolLeft,
unsigned*&leftIbegin):

x(x),
y(y),
NofLocation(NofLocation),
//kernel(kernel),
//weightEPS(weightEPS),
sd(sd),
sumVar(sumVar),
outputCorM(outputCorM),
corM(corM),
indexV(indexV),
//ind(ind),
//w(w),
//kernelConvCovar(kernelConvCovar),
//locationKernelSD(locationKernelSD),
localGkernel(localGkernel),
covarSum(covarSum),
NofColumnProcessed(NofColumnProcessed),
NofCore(NofCore),
someoneQuit(someoneQuit),
NcolLeft(NcolLeft),
leftIbegin(leftIbegin){}

//when writing operator() function, it must use std::size_t!! Anything else won't pass!!
void operator()(std::size_t st, std::size_t end){
for(std::size_t I=st;I!=end;++I)
{
unsigned avgColForEachCore=NcolLeft/NofCore;
unsigned indexVi=I*avgColForEachCore, indexViend=indexVi+avgColForEachCore;

if(I==NofCore-1)indexViend=NcolLeft;

//unsigned kernelSize=kernel.size();
valtype covar=0;

std::clock_t timestart=std::clock();

for(;;++indexVi,++NofColumnProcessed[I])
{
  if(someoneQuit||indexVi==indexViend)
  {
    leftIbegin[I]=indexVi;
    break;
  }

  unsigned i=indexV[indexVi];

  std::clock_t timeend=std::clock();
  if((timeend-timestart)/(double)CLOCKS_PER_SEC>300)
  {
    std::cout<<-std::accumulate(NofColumnProcessed, NofColumnProcessed+NofCore, 0);
    timestart=timeend;
  }

  valtype covar2=0;
  for(unsigned j=i+1,jend=NofLocation;j!=jend;++j)
  {
    valtype corr=kernelcorr(localGkernel[i], localGkernel[j]);
    covar2+=corr*(sd[i]*sd[j]);


    // if(std::isnan(covar2)||std::isinf(covar2))
    // {
    //   std::ofstream myfile;
    //   myfile.open("C:/Users/i56087/Desktop/validation/NShurricane/error.csv",
    //               std::ofstream::out|std::ofstream::app);
    //   myfile<<i<<","<<j<<std::endl;
    //   myfile.close();
    // }



    if(outputCorM)
    {
      corM[NofLocation*i+j]=corr;
      corM[NofLocation*j+i]=corr;
    }
  }
  covar+=covar2;
}

covarSum[I]+=covar;
someoneQuit=1;
}
}
};






valtype localGtotalVarParallel(double*x, double*y, unsigned NofLocation,
                 std::vector<G2>&kernel, double weightEPS,
                 double*sd, valtype sumVar, bool outputCorM, double*corM,
                 unsigned NofCore, unsigned*index,
                 std::vector<unsigned>&indexVec, bool useMahalanobisDforWeight){

std::vector<std::vector<unsigned> >ind;
std::vector<std::vector<valtype> >w;
locationKernelIndWeight(x, y, NofLocation, ind, w, kernel, weightEPS, useMahalanobisDforWeight);

//unsigned kernelSize=kernel.size();
std::vector<valtype>kernelConvCovar;

kernelcovar(kernel, kernelConvCovar);

std::vector<valtype>locationKernelSD;

computelocationKernelSD(x, y, NofLocation, kernel, weightEPS,
                        locationKernelSD, kernelConvCovar, ind, w);


// std::ofstream myfile;
// myfile.open("C:/Users/i56087/Desktop/validation/NShurricane/kernel.csv");

std::vector<G2>localG(NofLocation);
for(unsigned i=0;i!=NofLocation;++i)
{
  localG[i].formLocovar(x[i],y[i],&kernel.front(),ind[i],w[i]);
  // myfile<<localG[i].miu1<<","<<localG[i].miu2<<","<<localG[i].var1<<","
  //       <<localG[i].covar<<","<<localG[i].var2<<std::endl;
}

//myfile.close();

std::vector<valtype>covarVec(NofCore,0);
valtype*covarSum=&covarVec.front();

std::vector<unsigned>NofColumnProcessedV(NofCore,0);
unsigned*NofColumnProcessed=&NofColumnProcessedV.front();


{

  unsigned c1[NofCore], *coreStart=c1;
  unsigned c2[NofCore], *leftIbegin=c2;

  while(true)
  {
    unsigned*indexV=&indexVec.front();
    volatile bool someoneQuit=0;

    unsigned NcolLeft=indexVec.size();
    //if(NcolLeft==0)break;

    if(NcolLeft<NofCore)NofCore=1;

    unsigned avgLoadEachCore=NcolLeft/NofCore;
    for(unsigned i=0;i!=NofCore;++i)coreStart[i]=i*avgLoadEachCore;

    //std::cout<<NofCore<<" "<<NcolLeft<<"\n";

    NSredistributeLocovar T(x,y,NofLocation,sd,sumVar,
      outputCorM,corM,indexV,localG,covarSum,NofColumnProcessed,NofCore,
        someoneQuit,NcolLeft,leftIbegin);

    parallelFor(0, NofCore, T);

    //std::cout<<1.1<<std::endl;

    if(NofCore==1)break;

    for(int i=NofCore-1;i>=0;--i)
    {
      //std::cout<<"bounds=="<<coreStart[i]<<" "<<leftIbegin[i]<<std::endl;
      indexVec.erase(indexVec.begin()+coreStart[i], indexVec.begin()+leftIbegin[i]);
    }
    //std::cout<<1.2<<std::endl;

  }
}


valtype covar=accumulate(covarVec.begin(), covarVec.end(), 0.0);
if(outputCorM)
{
  for(unsigned i=0;i!=NofLocation;++i)corM[i*NofLocation+i]=1;
}
return sumVar+covar+covar;
}















// [[Rcpp::export]]
double GMKtotalVarCppParallel(List location, List GMkernel, NumericVector locationVar,
                      NumericMatrix corM, IntegerVector indexV, IntegerVector indexVector,
                      double weightEPS, bool outputCorM, //bool equalDivide,
                      bool absorbWeightInLocalCovar, bool useMahalanobisDforWeight){

std::vector<G2>kernel(GMkernel.size());
for(int i=0,iend=kernel.size();i!=iend;++i)
{
  NumericVector tmp=GMkernel[i];
  kernel[i].w=tmp[0];
  kernel[i].miu1=tmp[1];
  kernel[i].miu2=tmp[2];
  kernel[i].var1=tmp[3];
  kernel[i].covar=tmp[4];
  kernel[i].var2=tmp[5];
}
NumericVector x=location[0], y=location[1];

std::vector<valtype>locationSD(locationVar.size());
for(unsigned i=0,iend=locationSD.size();i!=iend;++i)
  locationSD[i]=std::sqrt(locationVar[i]);

unsigned NofCore=indexV.size();
unsigned index[NofCore];

std::copy(indexV.begin(),indexV.end(),index);

std::vector<unsigned>indexVec;
//if(!equalDivide)
indexVec.assign(indexVector.begin(), indexVector.end());

if(!absorbWeightInLocalCovar)return totalVarParallel(&x[0], &y[0], locationVar.size(),
   kernel, weightEPS, &locationSD.front(), std::accumulate(locationVar.begin(),
    locationVar.end(),0.0), outputCorM, &*corM.begin(),
      indexV.size(), index, indexVec, useMahalanobisDforWeight);


return localGtotalVarParallel(&x[0], &y[0], locationVar.size(),
        kernel, weightEPS, &locationSD.front(),
          std::accumulate(locationVar.begin(),locationVar.end(),0.0),
            outputCorM, &*corM.begin(), indexV.size(), index, indexVec, useMahalanobisDforWeight);


}

































