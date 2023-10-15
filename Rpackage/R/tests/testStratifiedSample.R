

Rcpp::sourceCpp('src/tests/testStratifiedSample.cpp', cacheDir = '../tempFiles',
                verbose = T)
repl = sample(c(0, 1), 1)
Npopu = 100000
sampleSize = 80000
sed = sample(1e6, 1)
system.time({
tmp = testStratifiedSampling(Npopu, sampleSize, replace = repl, seed = sed)
})
if (!repl && length(unique(tmp)) != sampleSize)
  cat("!repl && length(unique(tmp)) != sampleSize")


sed = 1:10000
tmp = lapply(sed, function(x) testStratifiedSampling(3, 2, replace = T, seed = x))
table(unlist(lapply(tmp, function(x) paste0(x, collapse = '-'))))


