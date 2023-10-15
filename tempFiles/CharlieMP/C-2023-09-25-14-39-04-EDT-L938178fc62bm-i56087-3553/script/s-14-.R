tmpDir = '/finance_develop/Charlie/BayesianFitDistr/tempFiles/CharlieMP/C-2023-09-25-14-39-04-EDT-L938178fc62bm-i56087-3553'
i = 14L
lg = file(paste0(tmpDir, '/log/t-', i, '-.txt'), open = 'wt')
sink(lg, type = 'output')
sink(lg, type = 'message')
load(paste0(tmpDir, '/input/t-', i, '-.Rdata'))
load(paste0(tmpDir, '/commonData.Rdata'))
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
