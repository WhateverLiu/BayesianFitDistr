\name{convBtRgd}
\alias{convBtRgd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Brute force convolve and then regrid
}
\description{
Convolve 2 distributions with brute force and regrid the result to user-defined grids.
}
\usage{
convBtRgd(X, Y, rgMethod="lr", N = 64L)
}
\arguments{
\item{X}{A 2-column data frame as the 1st distribution. 1st column is the variable values, 2nd column is the associated probabilities.
}
\item{Y}{
A 2-column data frame as the 2nd distribution.
}
\item{rgMethod}{
Regridding method. "lr" denotes linear regridding, "r" denotes new regridding
}
\item{N}{
Number of grids for regridding the resulting distribution. Default to 64.
}
}
\value{A 2-column data frame as the distribution of the convolution of X and Y.
}
\details{X and Y should have been sorted by their 1st column before going in.}

