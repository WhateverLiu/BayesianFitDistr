\name{armaSolve}
\alias{armaSolve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solve a linear system
}
\description{
A simple interface to Armadillo's dense solver.
}
\usage{
armaSolve(X, y, fastsolve = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Observation matrix.
}
  \item{y}{
Response values.
}
  \item{fastsolve}{
Switch on \code{arma::solve_opts::fast} if true.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Coefficient vector.
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
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
