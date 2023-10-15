\name{geocodeToGridIndex}
\alias{geocodeToGridIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate grid cell index
}
\description{
Gnerate grid cell index given latitude and longitude
}
\usage{
geocodeToGridIndex(geocode, hierachySize = c(1/120, 5/120), startLat = -90, startLong = -180)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geocode}{
A 2-column data frame. Every row is a claim's (longitude, latitude)
}
  \item{hierachySize}{
Latitude and longitude's delta for each hierachy. The latter hierachy's delta must be a multiple of its prior.
}
  \item{startLat}{
The latitude of the 1st grid cell's left lower corner. AIR set the point (lat= -90, long= -180) as the starting location for indexing the grid cells.
}
  \item{startLong}{
The longitude of the 1st grid cell's left lower corner.
}
}
\value{
The grid cell index. Number of columns is twice the length of "hierachySize".
}
\details{
As an example, consider a claim located at (latitude=30.38261, longitude=50.86317), then its first hierachy grid cell index is (int((30.38261-(-90))/(1/120))+1, int((50.86317-(-180))/(1/120))+1); its second hierachy grid cell index is (int((30.38261-(-90))/(5/120))+1, int((50.86317-(-180))/(5/120))+1)
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
To use gridIndex, one should first sort it by from the last column to the first column
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
