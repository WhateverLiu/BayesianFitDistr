\name{plotCoor}
\alias{plotCoor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Produce actual coordinates from relative coordinate in plot.
}
\description{
Given the coordinate in \code{[0, 1] x [0, 1]}, generate the actual coordinate for placing texts or legends in figure.
}
\usage{
plotCoor(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x coordinate of the relative location in plot. Should be in [0, 1].
}
  \item{y}{
y coordinate of the relative location in plot. Should be in [0, 1].
}
}
\details{
A \code{plot()} function must have been called.
}
\value{
A numeric vector of size 2: \code{(x, y)}
}
