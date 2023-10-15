\name{generateCorMultiplier}
\alias{generateCorMultiplier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate the multipliers on correlation coefficients.
}
\description{
The inner product of the result and the corresponding correlation coefficient vector will be the theoretical variance of the aggregated distribution.
}
\usage{
generateCorMultiplier(var, gridIndR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var}{
A numeric vector. Each element is the variance of a distribution.
}
  \item{gridIndR}{
Equivalent to "gridIndex" in mixConvSpatialClaimSampleList(). A data frame where each row is the grid index of the distribution corresponding to "var".
}
}
\details{
For example, if this function returns a vector (50, 100, 200, 300), and we want to impose correlation (0.06, 0.03, 0.01) on each grid level, then the variance of the aggregated distribution will be sum(c(50, 100, 200, 300)*c(1, 0.06, 0.03, 0.01)).

This function comes in handy for inverse engineering (optimizing) correlation coefficients.
}
\value{
A numeric vector as the multipliers. Its length would be equal to ncol(gridIndR)/2+1
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
