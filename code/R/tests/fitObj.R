

Rcpp::sourceCpp('cpp/amalgamation.cpp')


l = 1
u = 20
abc = sort(runif(3, l, u))
abc[1:3] = sort(abc[1:3])
abc[1:2] = sample(abc[1:2])
lm1 = runif(1)
abc_lm1 = c(abc, lm1)
d = solve_d(abc_lm1, 1e-8, 100, T)
abcd = c(abc, d)
actuarLevtrbeta(abcd); lm1
drs = 
drs = empc(actuarRtrbeta(1000, abcd), n = 64)




tmp = fitObj(x, p, size, abc_l1m, eps, maxit, method = "llh")























