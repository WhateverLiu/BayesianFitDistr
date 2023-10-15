\name{sourceCppExtra}
\alias{sourceCppExtra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Enhanced Rcpp::sourceCpp()
}
\description{
Enhanced Rcpp::sourceCpp(). The function can accept .o and .dll files to compile.
}
\usage{
sourceCppExtra(
  file = "", code = NULL,
  env = globalenv(), embeddedR = TRUE,
  rebuild = FALSE,
  cacheDir = getOption("rcpp.cache.dir", tempdir()),
  cleanupCacheDir = FALSE,
  showOutput = verbose,
  verbose = TRUE,
  dryRun = FALSE,
  windowsDebugDLL = FALSE,
  echo = TRUE,
  externObjLibFiles = NULL,
  retrieveSLL = FALSE,
  dllsNeeded = NULL
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
A character string giving the path name of a file
}
  \item{code}{
See the same argument in Rcpp::sourceCpp()
}
  \item{env}{
See the same argument in Rcpp::sourceCpp()
}
  \item{embeddedR}{
See the same argument in Rcpp::sourceCpp()
}
  \item{rebuild}{
See the same argument in Rcpp::sourceCpp()
}
  \item{cacheDir}{
See the same argument in Rcpp::sourceCpp()
}
  \item{cleanupCacheDir}{
See the same argument in Rcpp::sourceCpp()
}
  \item{showOutput}{
See the same argument in Rcpp::sourceCpp()
}
  \item{verbose}{
See the same argument in Rcpp::sourceCpp()
}
  \item{dryRun}{
See the same argument in Rcpp::sourceCpp()
}
  \item{windowsDebugDLL}{
See the same argument in Rcpp::sourceCpp()
}
  \item{echo}{
See the same argument in Rcpp::sourceCpp()
}
  \item{externObjLibFiles}{
A character vector of full paths of the .o or .lib files.
}
  \item{retrieveSLL}{
TRUE to copy back the compiled .o file into the directory where \code{file} resides.
}
  \item{dllsNeeded}{
A character vector of full paths of .dll files.
}
}
\details{
}
\value{
See the same argument in Rcpp::sourceCpp().
}

\examples{
# sourceCppExtra("externLib/test.cpp", rebuild = T,
#                verbose = T, retrieveSLL = TRUE)
# 
# # There should already have been a .h file
# # included in testLinkLib.cpp. This .h file
# # declares functions declared in test.cpp
# sourceCppExtra("testLinkLib.cpp", verbose = T,
#                rebuild = T, externObjLibFiles = "externLib/test.o")
# x = runif(10); y = runif(10)
# range(timesTwo(x, y) - (x + y) * x)
}



