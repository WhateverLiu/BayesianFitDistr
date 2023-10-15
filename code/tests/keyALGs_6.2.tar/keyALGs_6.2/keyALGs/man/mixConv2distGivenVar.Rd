\name{mixConv2distGivenVar}
\alias{mixConv2distGivenVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mix-convolve two distributions given the variance
}
\description{
Mix-convolve 2 distributions given the wanted variance of the mix-convolution
}
\usage{
mixConv2distGivenVar(dist1, dist2, varWanted, regridMethod = 'r', N = 256L, dist1var = -1, dist2var = -1, useFFT = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist1}{
A 2-column data frame as the 1st distribution
}
  \item{dist2}{
A 2-column data frame as the 2nd distribution
}
  \item{varWanted}{
The wanted variance
}
  \item{regridMethod}{
Regridding method fed to the convolution algorithm. "l"(first letter of "linear") denotes linear regridding. "r" denotes new regridding
}
  \item{N}{
Maximal number of grid points in the resultant distribution
}
  \item{dist1var}{
The variance of the 1st distribution. If not given (-1), the function will compute first.
}
  \item{dist2var}{
The variance of the 2nd distribution. If not given (-1), the function will compute first.
}
\item{useFFT}{
Whether use FFT to convolve.
}}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A 2-column data frame as the mix-convolution
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
