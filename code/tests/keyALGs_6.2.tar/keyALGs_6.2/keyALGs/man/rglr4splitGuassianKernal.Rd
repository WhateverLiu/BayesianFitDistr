\name{rglr4splitGuassianKernal}
\alias{rglr4splitGuassianKernal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Four dimensional regrid with Gaussian kernal
}
\description{
A modified four dimensional regridding method
}
\usage{
rglr4splitGuassianKernal(X, ngd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A 2-column data frame as the distribution
}
  \item{ngd}{
A numeric vector as the new grid
}
}
\details{
The 4 dimensional regridding method uses triangular kernal. This method uses Guassian kernal. It is for research purpose.
}
\value{
A 2-column data frame as the new distribution
}
