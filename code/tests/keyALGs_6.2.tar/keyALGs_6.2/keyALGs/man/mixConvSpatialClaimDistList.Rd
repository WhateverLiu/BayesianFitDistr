\name{mixConvSpatialClaimDistList}
\alias{mixConvSpatialClaimDistList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mix-convolve a list of distributions
}
\description{
Mix-convolve a list of distributions that are correlated given a multi-hierarchical spatial correlation model
}
\usage{
mixConvSpatialClaimDistList(distList, gridIndex, cors, regridMethod = "r", headTrunc = 1e-100, tailTrunc = 1e-100, N = 256L, forDeltaCompare = 256L, single = 1, useFFT = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distList}{
A list of distributions. Each distribution is a 2-column data frame. 1st column is variable values, 2nd column is probabilities
}
  \item{gridIndex}{
Please refer to "gridIndex" in mixConvSpatialClaimSampleList()
}
  \item{cors}{
Please refer to "corsR" in mixConvSpatialClaimSampleList()
}
  \item{regridMethod}{
The regridding method employed inside the algorithm. "l" (first letter of "linear") denotes linear regridding, "r" denotes the new regridding
}
  \item{headTrunc}{
Head truncation threshold fed to the convolution algorithm
}
  \item{tailTrunc}{
Tail truncation threshold fed to the convolution algorithm
}
  \item{N}{
Number of grid points in resulting distribution fed to the convolution algorithm
}
  \item{forDeltaCompare}{
Parameter for delta comparison fed to the convolution algorithm
}
  \item{single}{
Truncation method fed to the convolution algorithm
}
\item{useFFT}{
Whether to use FFT for convolution}
}
\details{
This function is just a wrapper of function "mixConvSpatialClaimSampleList" with arguments sampleListR=list(1:length(distList)), fullDistListR=distList, cors=corsR, maxCore=1
}
\value{
A list of size 2. 1st element is the theoretical variance supposed to have, 2nd element is the resultant mix-convolution
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
On the console the function will print the total number of fails to impose the supposed variance.

In distList, a NULL distribution, or a distribution having no variance, should have only 1 row. It is also recommended to remove all the NULL distributions prior to calling the function.
}

\examples{
N=20000
geoCode=data.frame(long=runif(N, -120, -70), lat=runif(N, 25, 48))


# make 20000 distributions for the locations:
# distribution loss min
Min=runif(N, 0, 1e6)
Max=runif(N, 0, 1e6)
Len=sample(2L:100L, N, replace = TRUE)
# don't accept distribution of length 1. Please
# make sure your distributions correct

distList=mapply(function(x,y,z)
{
  val=seq(min(x,y), max(x,y), len=z)
  P=abs(rnorm(z))
  P=P/sum(P)
  data.frame(val, P)
}, Min, Max, Len, SIMPLIFY = FALSE)


# provide 5 spatial level: 1km cell, 5km cell, 20km cell, 50km cell, global
gridInd=geocodeToGridIndex(geoCode, hierachySize = c(1/120, 5/120, 20/120, 50/120))
gridInd$long5=0L; gridInd$lat5=0L # the global scale


# order the locations and distributions from the largest scale to the smallest scale
tmpOrder=do.call(order, gridInd[ncol(gridInd):1L])
gridInd=gridInd[tmpOrder,]
distList=distList[tmpOrder]
geoCode=geoCode[tmpOrder,] # unnecessary


# provide the correlations on each scale:
Cor=c(0.5,0.4,0.3,0.2,0.1)


# mix-convolve the distributions:
# use 4-point regridding, N=256L(distribution length mustn't be greater
# than 256 at any time point including the initial distributions, as
# the software group does. Don't change it if unfamiliar with the
# concepts), tail truncation at 1e-10 (software group does), head
# truncation at 1e-10 (software group don't do head truncation which
# could hurt precision if linear regridding and naive sequential order of
# convolution is used. This function implements the best order of
# convolution/recursive nearest neighbor so truncation and regridding
# method wouldn't matter that much. If truncation not wanted, set
# headTrunc=1e-100). The rest parameters are for research purpose and
# should only be changed by expert.

aggregateDist=mixConvSpatialClaimDistList(distList, gridIndex=gridInd, cors=Cor, regridMethod="r", N=256L, tailTrunc=1e-10, headTrunc=1e-100, forDeltaCompare=256L, single = 1L, useFFT=0)

# check if the theoretical variance is imposed
Var(aggregateDist[[2]])/aggregateDist[[1]]-1

# plot the distribution
plot(aggregateDist[[2]],type="h")

# Another version of CLT: convolve finite discrete distributions (don't
# have to be identical) will converge to Gaussian. So mix-convolving them
# with small correlations will still result in Gaussian-shape
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
