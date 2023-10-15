\name{covarMatSum}
\alias{covarMatSum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the theoretical variances for subportfolios
}
\description{
Given a multi-hierarchical spatial correlation model, compute the variances supposed to be imposed for a list of subportfolios in a parallelizable fashion
}
\usage{
covarMatSum(sampleList, variances, cors, gridInd, maxCore = 8L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleList}{
A list of integer vectors. Each integer vector MUST BE SORTED. Each integer vector, as a vector of indices, denotes a subportfolio from "variances". Please refer to the same argument in function "mixConvSpatialClaimSampleList" for more information
}
  \item{variances}{
A numeric vector of variances belonging to each distribution in the the entire distribution portfolio.
}
  \item{cors}{
The same argument as corsR in mixConvSpatialClaimSampleList()}
  \item{gridInd}{
The same argument as gridIndex in mixConvSpatialClaimSampleList()
}
  \item{maxCore}{
Number of CPUs wanted for parallel computing. MUST NOT BE GREATER THAN THE ACTUAL NUMBER OF IDLE CPUs ON MACHINE. MUST NOT BE GREATER THAN THE LENGTH OF "variances"
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A numeric vector. The i_th element is the theoretical variance for the i_th subportfolio.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
If parametric loss aggregation is of interest, this function is used to calculate the aggregated distribution's variance.
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
