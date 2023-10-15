\name{paraOLS}
\alias{paraOLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parallelized least square fitting for dense observation matrix
}
\description{
The speedup comes from multithreading the Gram matrix computation.
}
\usage{
paraOLS(
  y, X,
  maxCore = 7L,
  verboseInt = 8192L,
  useGaussSeidel = FALSE,
  maxIter = 1000000L,
  convergenceEPS = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Observed response vector.
}
  \item{X}{
Observation matrix.
}
  \item{maxCore}{
Maximal threads to invoke.
}
  \item{verboseInt}{
Print progress every \code{verboseInt} inner products.
}
  \item{useGaussSeidel}{
Use Gauss-Seidel iteration for solving the Gramian matrix. If \code{X}'s columns are linearly independent, then \code{X}'\code{X} is positive definite, which guarantees Gauss-Seidel's convergence.
}
  \item{maxIter}{
Maximal iterations for Gauss-Seidel.
}
  \item{convergenceEPS}{
Convergence criterion for Gauss-Seidel.
}
}

\value{

A list of 3

  coef    : num [1:101, 1] 0.0619 0.0526 ...   = Coefficients.

  pred    : num [1:1000] -0.04374 0.43178 ...  = Model responses.

  totalErr: num 902 = Sum of squared error.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
