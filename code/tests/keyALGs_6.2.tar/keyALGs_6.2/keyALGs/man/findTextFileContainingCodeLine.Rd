\name{findTextFileContainingCodeLine}
\alias{findTextFileContainingCodeLine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the text code file that contains a certain code line
}
\description{
Find the text code file that contains a certain code line.
}
\usage{
findTextFileContainingCodeLine(
  direct,
  codeLine,
  filetype = c("R", "r"),
  findAll = F,
  recursive = T,
  readCharNchar = 1e7L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{direct}{
Directory to search in.
}
  \item{codeLine}{
The code line.
}
  \item{filetype}{
Target file extension.
}
  \item{findAll}{
TRUE to find all files that contain the code line. FALSE to find the first
file that contains the code line.
}
  \item{recursive}{
TRUE to search in all subdirectories.
}
\item{readCharNchar}{
Number of characters to read when using readChar().
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of file names.
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
