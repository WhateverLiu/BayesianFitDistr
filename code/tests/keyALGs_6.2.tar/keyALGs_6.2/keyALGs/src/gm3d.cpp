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

#define vec std::vector
#define valtype double
#define indtype unsigned

valtype sqrt8PI3=std::sqrt(8*M_PI*M_PI*M_PI);


struct dynamicTasking{
unsigned NofCore;
unsigned NofAtom;
unsigned*current,*UB;
std::vector<unsigned>container;
tbb::mutex m;

dynamicTasking(){}

void reset(unsigned NofCPU, unsigned NofTask, bool atomize){
NofCore=NofCPU;
NofAtom=NofTask;
}

dynamicTasking(unsigned NofCPU, unsigned NofTask, bool atomize){
reset(NofCPU, NofTask, atomize);
}

bool nextTaskID(unsigned&taskID)
{
  bool rst=1;
  m.lock();
  if(NofAtom==0)rst=0;
  else
  {
    taskID=NofAtom-1;
    --NofAtom;
  }
  m.unlock();
  return rst;
}

void reset(unsigned NofCPU, unsigned NofTask)
{
  container.resize(2*NofCPU);
  current=&*container.begin();
  UB=&*container.begin()+NofCPU;
  unsigned avgTask=NofTask/NofCPU;
  current[0]=0;
  UB[0]=current[0]+(NofTask-NofCPU*avgTask)+avgTask;
  for(unsigned i=1,iend=NofCPU;i!=iend;++i)
  {
    current[i]=UB[i-1];
    UB[i]=current[i]+avgTask;
  }
  NofCore=NofCPU;
}

dynamicTasking(unsigned NofCPU, unsigned NofTask)
{
  reset(NofCPU, NofTask);
}

bool nextTaskID(unsigned cpuI, unsigned&taskID)
{
  //what u gonna do once this chunk is finished
  if(current[cpuI]==UB[cpuI])
  {
    m.lock();
    unsigned maxleft=0, whichone=0;
    for(unsigned i=0;i!=NofCore;++i)
    {
      if(maxleft<UB[cpuI]-current[cpuI])
      {
        maxleft=UB[cpuI]-current[cpuI];
        whichone=i;
      }
    }
    if(maxleft<=1)
    {
      m.unlock();
      return 0;
    }
    UB[cpuI]=UB[whichone];
    current[cpuI]=current[whichone]+(UB[whichone]-current[whichone])/2;
    UB[whichone]=current[cpuI];
    taskID=current[cpuI];
    ++current[cpuI];
    m.unlock();
  }
  else
  {
    taskID=current[cpuI];
    ++current[cpuI];
  }
  return 1;
}
};


struct g3d{
valtype mu0, mu1, mu2, C00, C01, C02, C11, C12, C22, alpha;
vec<valtype>ptr;

valtype eval(valtype x0, valtype x1, valtype x2, bool withAlpha=1)
{
  valtype A=C11*C22-C12*C12, B=C01*C22-C02*C12, C=C01*C12-C02*C11,
    det=C00*A-C01*B+C02*C,
    z0=x0-mu0, z1=x1-mu1, z2=x2-mu2;

  valtype tmp=A*z0*z0+(C00*C22-C02*C02)*z1*z1+(C00*C11-C01*C01)*z2*z2-
    2*(B*z0*z1-C*z0*z2+(C00*C12-C01*C02)*z1*z2);

  tmp=std::exp(-tmp/(det*2))/(sqrt8PI3*std::sqrt(det));
  if(std::isnan(tmp))return 0;
  if(withAlpha)return tmp*alpha;
  return tmp;
}
};


struct val3d{valtype x0, x1, x2;};


// [[Rcpp::export]]
double testG3dVal(NumericVector x, NumericVector para)
{
  g3d G;
  valtype*Gptr=(valtype*)(&G);
  std::copy(&*para.begin(),&*para.end(),Gptr);
  return G.eval(x[0],x[1],x[2]);
}


struct cmptDensity: public Worker{
g3d*gmodel;
indtype gmodelSize;
val3d*X;
indtype Xsize;
dynamicTasking dT;


void computeDensity(unsigned i, g3d*gmodel, val3d*X, indtype Xsize)
{
  // iend is Xsize*gmSize
  indtype whichModel=(indtype)i/Xsize, offset=i % Xsize;
  g3d&gm=gmodel[whichModel];
  gm.ptr[offset]=gm.eval(X[offset].x0, X[offset].x1, X[offset].x2);
}


void operator()(std::size_t st, std::size_t end){
for(std::size_t I=st;I!=end;++I)
{
  for(;;)
  {
    unsigned objI=0;
    if(!dT.nextTaskID(I, objI))break;
    computeDensity(objI, gmodel, X, Xsize);
  }
}
}

cmptDensity(g3d*gmodel, indtype gmodelSize, val3d*X, indtype Xsize, unsigned NofCPU):
gmodel(gmodel),gmodelSize(gmodelSize),X(X),Xsize(Xsize)
{
  dT.reset(NofCPU, gmodelSize*Xsize);
  parallelFor(0, dT.NofCore, *this);
}
};


struct cmptRowSum: public Worker{
valtype*rst;
g3d*gmodel;
indtype gmodelSize;
indtype Xsize;
dynamicTasking dT;


void computeRowSum(valtype*rst, indtype rowi, g3d*gmodel, indtype gmodelSize)
{
  rst[rowi]=0;
  for(indtype i=0,&iend=gmodelSize;i!=iend;++i)
  {
    rst[rowi]+=gmodel[i].ptr[rowi];
  }
}


void operator()(std::size_t st, std::size_t end){
for(std::size_t I=st;I!=end;++I)
{
  for(;;)
  {
    unsigned objI=0;
    if(!dT.nextTaskID(I, objI))break;
    computeRowSum(rst, objI, gmodel, gmodelSize);
  }
}
}

cmptRowSum(valtype*rst, g3d*gmodel, indtype gmodelSize, indtype Xsize, unsigned NofCPU):
gmodel(gmodel),gmodelSize(gmodelSize),Xsize(Xsize)
{
  dT.reset(NofCPU, Xsize);
  parallelFor(0, dT.NofCore, *this);
}
};


struct cmptWeight: public Worker{
g3d*gmodel;
indtype gmodelSize;
val3d*X;
indtype Xsize;
valtype*rowSumV;
dynamicTasking dT;


void singleWeight(indtype i, valtype*rowSumV, g3d*gmodel, indtype gmodelSize, indtype Xsize)
{
  indtype whichModel=(indtype)i/Xsize, offset=i%Xsize;
  gmodel[whichModel].ptr[offset]/=rowSumV[offset];
}


void operator()(std::size_t st, std::size_t end){
for(std::size_t I=st;I!=end;++I)
{
  for(;;)
  {
    unsigned objI=0;
    if(!dT.nextTaskID(I, objI))break;
    singleWeight(objI, rowSumV, gmodel, gmodelSize, Xsize);
  }
}
}

cmptWeight(g3d*gmodel, indtype gmodelSize, val3d*X, indtype Xsize, unsigned NofCPU):
  gmodel(gmodel), gmodelSize(gmodelSize), X(X), Xsize(Xsize)
{
  dT.reset(NofCPU, gmodelSize*Xsize);
  parallelFor(0, dT.NofCore, *this);
}
};


struct updatePara: public Worker{
g3d*gmodel;
indtype gmSize;
val3d*X;
indtype Xsize;
dynamicTasking dT;

valtype update1gm(g3d&gm, val3d*X, indtype Xsize)
{
  valtype paraResv[10];
  std::copy((valtype*)(&gm), (valtype*)(&gm)+10, paraResv);

  //mu0, mu1, mu2, C00, C01, C02, C11, C12, C22, alpha;


  valtype Nk=0;
  gm.mu0=0;
  gm.mu1=0;
  gm.mu2=0;
  for(indtype i=0,&iend=Xsize;i!=iend;++i)
  {
    Nk+=gm.ptr[i];
    gm.mu0+=gm.ptr[i]*X[i].x0;
    gm.mu1+=gm.ptr[i]*X[i].x1;
    gm.mu2+=gm.ptr[i]*X[i].x2;
  }
  gm.alpha=Nk/Xsize;
  gm.mu0/=Nk;
  gm.mu1/=Nk;
  gm.mu2/=Nk;


  gm.C00=0;gm.C01=0;gm.C02=0;gm.C11=0;gm.C12=0;gm.C22=0;
  for(indtype i=0,&iend=Xsize;i!=iend;++i)
  {
    valtype x0_mu0=X[i].x0-gm.mu0, x1_mu1=X[i].x1-gm.mu1, x2_mu2=X[i].x2-gm.mu2;
    gm.C00+=x0_mu0*x0_mu0*gm.ptr[i];
    gm.C01+=x0_mu0*x1_mu1*gm.ptr[i];
    gm.C02+=x0_mu0*x2_mu2*gm.ptr[i];
    gm.C11+=x1_mu1*x1_mu1*gm.ptr[i];
    gm.C12+=x1_mu1*x2_mu2*gm.ptr[i];
    gm.C22+=x2_mu2*x2_mu2*gm.ptr[i];
  }
  gm.C00/=Nk; gm.C01/=Nk; gm.C02/=Nk; gm.C11/=Nk; gm.C12/=Nk; gm.C22/=Nk;


  valtype paraDiff=0, *tmpP=(valtype*)(&gm);
  for(indtype i=0,iend=10;i!=iend;++i)
  {
    valtype tmp=std::abs((paraResv[i]-tmpP[i])/(paraResv[i]+tmpP[i]));
    if(std::isnan(tmp))tmp=0;
    if(paraDiff<tmp)paraDiff=tmp;
  }
  return paraDiff;
}


void operator()(std::size_t st, std::size_t end){
for(std::size_t I=st;I!=end;++I)
{
  for(;;)
  {
    unsigned objI=0;
    if(!dT.nextTaskID(I, objI))break;
    update1gm(gmodel[objI], X, Xsize);
  }
}
}


updatePara(g3d*gmodel, indtype gmSize, val3d*X, indtype Xsize, unsigned NofCPU):
  gmodel(gmodel), gmSize(gmSize), X(X), Xsize(Xsize)
{
  dT.reset(NofCPU, gmSize);
  parallelFor(0, dT.NofCore, *this);
}
};






















































