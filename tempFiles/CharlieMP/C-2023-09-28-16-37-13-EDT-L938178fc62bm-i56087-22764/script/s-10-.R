tmpDir = '/finance_develop/Charlie/BayesianFitDistr/tempFiles/CharlieMP/C-2023-09-28-16-37-13-EDT-L938178fc62bm-i56087-22764'
i = 10L
lg = file(paste0(tmpDir, '/log/t-', i, '-.txt'), open = 'wt')
sink(lg, type = 'output')
sink(lg, type = 'message')
load(paste0(tmpDir, '/input/t-', i, '-.Rdata'))
load(paste0(tmpDir, '/commonData.Rdata'))
source('/finance_develop/Charlie/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_36cffa84562/P0ensembleUsingLinearInterpo.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c71d34eb24/assignP0.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c71962cadb/autoCorr.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c7113195b1a/correctMeanBias.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c71af33526/distanceAPI.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c7156961cae/extractMain.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c715d8c81e9/fineTuneDiscretizationFulfillQArequirements-002.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c7114a3691b/fitObj.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c712ce13731/ignoranceScore.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c71263d19a1/interpolateEmpDistr.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c711a5d5df3/longestSubseq.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c7110d66d1c/makeEmpDistr.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c714919e179/mixDistr.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c71fb4c52d/movingAverage.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_a12f4c0002b9/P0ensembleUsingLinearInterpo.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c717482a152/solve-d.cpp.R')
source('/finance_develop/Charlie/BayesianFitDistr/tempFiles/sourceCpp-x86_64-redhat-linux-gnu-1.0.11/sourcecpp_4c714811e928/upscaleAndBoundPMF.cpp.R')
for (x in loadedLibNames) eval(parse(text = paste0('library(', x, ')')))
verbose = 0.01
rst = list()
if (verbose > 0.5 || verbose <= 0) rst = lapply(X, function(x) fun(x, commonData)) else {

      gap = max(1L, as.integer(round(length(X) * verbose)))
      cat('N(items) =', length(X), ': ')
      for (k in 1:length(X))
      {
        rst[[k]] = fun(X[[k]], commonData)
        if (k %% gap == 0L) cat(k, '')
      }
    }
save(rst, file = paste0(tmpDir, '/output/rst-', i, '-.Rdata'))
finishName = paste0('f-', Sys.info()['nodename'], '-pid', Sys.getpid(), '-', i)
file.create(paste0(tmpDir, '/complete/', finishName))
sink(type = 'message'); sink(type = 'output'); close(lg)
