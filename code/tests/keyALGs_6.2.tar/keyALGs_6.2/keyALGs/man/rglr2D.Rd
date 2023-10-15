\name{rglr2D}
\alias{rglr2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear regrid a 2D discrete distribution
}
\description{
Linear regrid 2 2D discrete distribution
}
\usage{
rglr2D(dist, xnew, ynew)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist}{
A 3-column data frame. 1st and 2nd columns are the X and Y coordinates of the variable value points. 3rd column is probabilities.
}
  \item{xnew}{
The new grid for X. Intervals within must be equal.
}
  \item{ynew}{
The new grid for Y. Intervals within must be equal.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A 3-column data frame. 1st and 2nd column are the X and Y coordinates of the variable value points. 3rd column is probabilities.
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
d=data.frame(x=rnorm(100),y=runif(100,-1,2),P=rnorm(100)^2)
d$P=d$P/sum(d$P)
# a 2D distribution is generated

dnew=rglr2D(d,seq(min(d$x),max(d$x),len=15),seq(min(d$y),max(d$y),len=10))
if (require("rgl"))
{
  plot3d(dnew,type="h",zlim=c(0,max(dnew[[3]])))
}

sum(dnew$P)
# check sum of probabilities

sum(dnew$P*dnew$x)-sum(d$P*d$x)
# check mean of x

sum(dnew$P*dnew$y)-sum(d$P*d$y)
# cehck mean of y
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line















