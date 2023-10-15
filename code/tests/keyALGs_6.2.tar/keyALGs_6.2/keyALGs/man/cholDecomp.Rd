\name{cholDecomp}
\alias{cholDecomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multithreaded Cholesky decomposition.
}
\description{
Multithreaded Cholesky decomposition. The source code provides a good example of how to remove the overhead of creating threads when parallelFor is inside a loop.
}
\usage{
cholDecomp(A, maxCore = 15L, useFullPara = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
A positive definite matrix.
}
  \item{maxCore}{
Maximum number of threads available.
}
  \item{useFullPara}{
FALSE activates the parallel version where parallelFor resides inside a loop. TRUE and FALSE comparison gives an idea of how much system overhead is incurred for forking and joining threads.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The upper triangle matrix from Cholesky decomposition.
}

