\name{normalizeDistributionWithAtoms}
\alias{normalizeDistributionWithAtoms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalize distributions such that its 1st moment matches the given means
}
\description{
This function will normalize the distribution by scaling. The method follows the section 4.4 post-optimization in the documentation of Fitting Fire Following Distribution.
}
\usage{
normalizeDistributionWithAtoms(X, support=NULL, min = 0, max, mean)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A list or a data frame of probabilities. If data frame, each column is a probability vector.
}
\item{support}{
A list or a data frame of distributions' supports. If data frame, each column is a distribution's support associated to X.
}
  \item{min}{
Minimum of the support if support is NULL
}
  \item{max}{
Maximum of the support if support is NULL
}
  \item{mean}{
Distribution mean to impose
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list or a data frame of the normalized probabilities
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
