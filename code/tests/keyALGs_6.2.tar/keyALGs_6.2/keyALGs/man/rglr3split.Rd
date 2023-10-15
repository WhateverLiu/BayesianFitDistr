\name{rglr3split}
\alias{rglr3split}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Three dimensional regrid
}
\description{
A new regridding method that can preserve the 0th, 1st and 2nd moment
}
\usage{
rglr3split(X, ngd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A 2-column data frame as the distribution
}
  \item{ngd}{
A numeric vector as the new grid}
}
\details{
This method is an introduction of the 4 dimensional regridding method. More details can be found in the research document.
}
\value{
A 2-column data frame as the new distribution.
}
