\name{inverseRglr}
\alias{inverseRglr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inverse linear regrid
}
\description{
Given a distribution, linear interpolate the probabilities onto a denser grid
}
\usage{
inverseRglr(X, ngd, normalize = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A 2-column data frame as the distribution; 1st column is variable values
}
  \item{ngd}{
A numeric vector as the denser grid. Should not be wider than X's support
}
  \item{normalize}{
Whether to normalize the regridded distribution (preserve the 0th moment)
}
}
\details{
Usually used in plotting X while the support requires more points. 1st moment of X won't be preserved
}
\value{
The regridded distribution as a 2-column data frame
}
\references{
%% ~put references to the literature/web site here ~
}
