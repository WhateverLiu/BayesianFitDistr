\name{win2lnxZipCopy}
\alias{win2lnxZipCopy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Zip and copy files from Windows to Linux cluster
}
\description{
Zip files, subdirectories or both in the current directory, then copy the zipped file to Linux cluster using Powershell, and then unzip the file in Linux environment. The goal is to accelerate files transmission between the local Windows machine and the remote Linux cluster.
}
\usage{
win2lnxZipCopy(
  winFiles,
  lnxDir,
  winLnxPathAdd = "//RSGrid",
  linuxGrid = "rscgrid1finance1",
  passwdPath = "data/passwd.Rdata",
  removeZipOnWin = TRUE,
  removeZipOnLnx = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{winFiles}{
A vector of names of files or subdirectories of interest on the local Windows.
}
  \item{lnxDir}{
Path of the directory on the Linux cluster.
}
  \item{winLnxPathAdd}{
The extra characters attached to paths of files on the Linux cluster.
}
  \item{linuxGrid}{
Linux cluster grid engine name.
}
  \item{passwdPath}{
Password .Rdata. Content should a string.
}
  \item{removeZipOnWin}{
After the copy is done, should the zip file on Windows be deleted?
}
  \item{removeZipOnLnx}{
After the copy is done, should the zip file on Linux be deleted?
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
None
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
# win2lnxZipCopy(winFiles = c("alldata.Rdata", "scripts"), lnxDir = "//RSGrid/finance_develop/Charlie/SVChristian005/data")
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
