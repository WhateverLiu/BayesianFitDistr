\name{gm2dEMparallel}
\alias{gm2dEMparallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multithreaded 2D Gaussian Mixture Modeling using adaptive Expectation Maximization Algorithm
}
\description{
Multithreaded 2D Gaussian Mixture Modeling using adaptive Expectation Maximization Algorithm
}
\usage{
gm2dEMparallel(xCoor, xval, gmWeight, gmMean, gmVar, weightEPS = 1e-20, 
  convergeEPS = 1e-06, maxit = 500L, showProgress = 5L, NofCore = 7L,
  useRelativeDiff = 1, ellipseAxisRatioThreshold = 14.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xCoor}{
A 2-column data frame. 1st column as x-coordinate, 2nd column as y-column
}
  \item{xval}{
Value at each location (coordinate). A numeric vector
}
  \item{gmWeight}{
A numeric vector. Initilization of the weights for the Gaussian components
}
  \item{gmMean}{
A 2-column data frame, containing the same number of rows as the length of gmWeight. Initilization of the means of the Gaussian components
}
  \item{gmVar}{
A list having the same size of the length of gmWeight. Each list is a numeric vector of 3 elements --- var1, covar, var2
}
  \item{weightEPS}{
Threshold of the Gaussian component weights. Any Gaussian component having weight less than weightEPS will be eliminated
}
  \item{convergeEPS}{
Convergence threshold. If the distance between the current parameters and the previous ones is less than convergeEPS, program stops.
}
  \item{maxit}{
Maximal number of iterations
}
  \item{showProgress}{
If showProgress = 5, console will display the total number of iterations and the converging difference each 5 iterations.
}
  \item{NofCore}{
Available CPU cores put in use
}
  \item{useRelativeDiff}{
Boolean value. Use relative difference for examining convergence
}
\item{ellipseAxisRatioThreshold}{
Threshold of the covariance ellipse's aspect ratio.
}
}
\details{
Not recommended for small dataset
}
\value{
A list of numeric vectors. Each vector is of size 6. 1st element is weight, 2nd - 3rd elements is mean of the Gaussian. 4th element is var1, 5th element is covariance, 6th element is var2
}
\references{
\href{http://dataclustering.cse.msu.edu/papers/mixtureTPAMI.pdf}{Unsupervised learning of finite mixture models}.
}
