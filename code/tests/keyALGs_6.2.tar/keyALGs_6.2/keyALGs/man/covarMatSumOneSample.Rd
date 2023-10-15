\name{covarMatSumOneSample}
\alias{covarMatSumOneSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the theoretical variance
}
\description{
Given a multi-hierarchical spatial correlation model, compute the variances supposed to be imposed for a portfolio of distributions
}
\usage{
covarMatSumOneSample(variances, gridInd, cors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variances}{
A numeric vector. The i_th element is the variance of the i_th distribution in the portfolio
}
  \item{gridInd}{
Please refer to the same argument in function "covarMatSum"
}
  \item{cors}{
Please refer to the same argument in function "covarMatSum"}
}
\details{
This function is just a wrapper of function "covarMatSum" with argument sampleList=1:length(variances) and maxCore=1
}
\value{
The theoretical variance
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
