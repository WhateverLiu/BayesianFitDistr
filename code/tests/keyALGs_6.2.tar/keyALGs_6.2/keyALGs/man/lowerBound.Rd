\name{lowerBound}
\alias{lowerBound}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binary search for the lower bounds
}
\description{
Given a vector \code{x} and a sorted vector \code{v}, find the 1-based index of the first element in \code{v} no less than \code{x[i]}.
}
\usage{
lowerBound(x, v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric vector.
}
  \item{v}{
The sorted vector.
}
}
\details{
This is just a wrapper of STL function lower_bound().
}
\value{
An integer vector as the searched values' indices, 1-based.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
If search fails, return the length of v plus 1
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
