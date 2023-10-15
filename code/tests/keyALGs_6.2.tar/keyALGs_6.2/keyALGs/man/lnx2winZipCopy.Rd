\name{lnx2winZipCopy}
\alias{lnx2winZipCopy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The counterpart of \code{win2lnxZipCopy}.
}
\description{
The counterpart of \code{win2lnxZipCopy}.
}
\usage{
lnx2winZipCopy(
  lnxFiles,
  winDir,
  winLnxPathAdd = "//RSGrid",
  linuxGrid = "rscgrid1finance1",
  passwdPath = "data/passwd.Rdata",
  removeZipOnWin = TRUE,
  removeZipOnLnx = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lnxFiles}{
A vector of names of files or subdirectories of interest on the Linux cluster.
}
  \item{winDir}{
Path of directory on the local Windows machine.
}
  \item{winLnxPathAdd}{
The extra characters attached to paths of files on the Linux cluster.
}
  \item{linuxGrid}{
See the same argument in \code{win2lnxZipCopy}.
}
  \item{passwdPath}{
See the same argument in \code{win2lnxZipCopy}.
}
  \item{removeZipOnWin}{
See the same argument in \code{win2lnxZipCopy}.
}
  \item{removeZipOnLnx}{
See the same argument in \code{win2lnxZipCopy}.
}
}
\examples{
# Copy files on Linux cluster to the current working directory on the local Windows machine.
# win2lnxZipCopy(lnxFiles = "//RSGrid/finance_develop/Charlie/SVChristian005/result", winDir = ".")
}











