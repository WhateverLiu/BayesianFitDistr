// [[Rcpp::depends(RcppParallel)]]
# include <RcppParallel.h>
# include <atomic>


struct dynamicTasking
{
  std::size_t NofCore;
  std::size_t NofAtom;
  std::atomic<std::size_t> counter;


  void reset(std::size_t &NofCPU, std::size_t NofTask)
  {
    NofCore = NofCPU > NofTask ? NofTask : NofCPU;
    NofCPU = NofCore;
    NofAtom = NofTask;
    counter = 0;
  }


  dynamicTasking(std::size_t NofCPU, std::size_t NofTask)
  {
    reset(NofCPU, NofTask);
  }


  bool nextTaskID(std::size_t &taskID, std::size_t increment = 1)
  {
    taskID = counter.fetch_add(increment);
    return taskID < NofAtom;
  }
};




// T is the job object.
// T must have the following member functions:
// 1. Njobs(void): return an integer, the total number of jobs.
// 2. void run(std::size_t i): run the i_th job.
// If jobs are heterogeneous, it would be better that lower IDs should take
// longer time to finish.
template<typename T>
struct ParaFor: public RcppParallel::Worker
{
  std::size_t grainSize, Njobs;
  dynamicTasking *dT;
  T *obj;
  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t I = 0;
      if(!dT->nextTaskID(I, grainSize)) break;
      for(std::size_t Iend = std::min(I + grainSize, dT->NofAtom); I < Iend; ++I)
        obj->run(I);
    }
  }
  ParaFor(T &X, std::size_t maxCore = 15, std::size_t grainSize = 1):
    grainSize(grainSize)
  {
    Njobs = X.Njobs();
    obj = &X;
    dynamicTasking dt(maxCore, Njobs); dT = &dt;
    parallelFor(0, maxCore, *this);
  }
};


































