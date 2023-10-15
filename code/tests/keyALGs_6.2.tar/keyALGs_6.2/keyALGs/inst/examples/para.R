

# One motivation behind para() is to allow multiprocessing over Rcpp functions
# compiled and linked in the current R session, without the necessity of
# tracking and recompilation of those functions in every other process.
cppFile = "
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double stu(NumericVector x) { return std::accumulate(x.begin(), x.end(), 0.0); }
"
writeLines(cppFile, con = "tmpCpp.cpp")
Rcpp::sourceCpp("tmpCpp.cpp")


# ==============================================================================
# Example I:
# ==============================================================================
X = lapply(1:10000, function(x) runif(1000))
y = runif(1000)
f = function(x, y) { stu(x) + sum(x * y) }


rstTruth = lapply(X, function(x) f(x, y))


rstPara = para(X, commonData = y, fun = f, maxNprocess = 100, wait = TRUE, 
               cleanTempFirst = TRUE, keepTempFolder = TRUE,
               RscriptExePath = NULL)


# Compare results:
cat("Max difference =", max(abs(unlist(rstTruth) - unlist(rstPara))), "\n")




# ==============================================================================
# Example II:
# ==============================================================================
g = function(x, y) { stu(x) } # Function should always have 2 arguments.
rstParaNoWait = para(X, commonData = NULL, fun = g, maxNprocess = 30, 
                     wait = FALSE, cleanTempFirst = TRUE, 
                     keepTempFolder = TRUE, RscriptExePath = NULL)
cat("Files that will have the results:\n\n"); for(u in rstParaNoWait) cat(u, "\n")











































