\name{convSplitAtom4productsAtomUntouched}
\alias{convSplitAtom4productsAtomUntouched}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Split-atom 4 products atom untouched
}
\description{
A fast convolution algorithm
}
\usage{
convSplitAtom4productsAtomUntouched(X, Y, regridMethod = "lr", N = 256L, forDeltaCompare = 256L, leftTrunc = 0, rightTrunc = 0, useFFT = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A 2-column data frame as the first distribution
}
  \item{Y}{
A 2-column data frame as the second distribution
}
  \item{regridMethod}{
A string indicating the regridding method in process.

"lr": linear regrid

"r3": 3 dimensional regrid

"r4": 4 dimensional regrid

"r4g": 4 dimensional regrid with Gaussian kernal
}
  \item{N}{
Maximum grid points allowed in the resulting distribution
}
  \item{forDeltaCompare}{
For computing the delta of the resulting convolution. Should not be greater than N.
}
  \item{leftTrunc}{
Tail truncation threshold
}
  \item{rightTrunc}{
Head truncation threshold
}
\item{useFFT}{
Whether use FFT to convolve}
}
\details{
Please refer to the instruction of convolution algorithm for more details
}
\value{
A 2-column data frame as the convolution distribution
}
