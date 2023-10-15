\name{findPlotCoor}
\alias{findPlotCoor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find real coordinates in plot
}
\description{
Given relative coordinates in [0, 1] x [0, 1], find real coordinates for plotting.
}
\usage{
findPlotCoor(relaX = 0.5, relaY = 0.5, realX = NULL, realY = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{relaX}{
x coordinate in [0, 1] x [0, 1].
}
  \item{relaY}{
y coordinate in [0, 1] x [0, 1].
}
  \item{realX}{
Real x coordinate. If not NULL, then realX and realY are used.
}
  \item{realY}{
Real y coordinate.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A numeric vector of size 2, (x, y).
}
