\name{snowCplyWindows}
\alias{snowCplyWindows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiprocessing using the snow package
}
\description{
Multiprocess computation using the snow package.
}
\usage{
snowCplyWindows(X, FUN, ..., maxCore = 7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A vector or a list of data.
}
  \item{FUN}{
Function to be applied to each element in \code{X}.
}
  \item{\dots}{
Extra argument to \code{FUN}.
}
  \item{maxCore}{
Maximal number of cores to be used.
}
}
\details{
The function first divides X equally into \code{maxCore} pieces using function \code{splitList()}, then calls functions in package \code{snow}.
}
\value{
A list of the same size of \code{X}.
}
