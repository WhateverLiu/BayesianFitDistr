\name{matmul}
\alias{matmul}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast matrix multiplication
}
\description{
Matrix multiplication faster than base and Armadillo under both single-thread and multithread setting.
}
\usage{
matmul(X, Y, maxCore = 16L, bidir = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A numeric matrix.
}
  \item{Y}{
A numeric matrix. Number of rows should equal X's number of columns.
}
  \item{maxCore}{
Number of threads available.
}
  \item{bidir}{
Bi-directional looping. An experiment for testing cache. \code{bidir = FALSE} should almost always be preferred.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A numeric matrix.
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
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
