\name{partition3wayLRMconv}
\alias{partition3wayLRMconv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Three-way partitioning accumulation
}
\description{
The order of the convolution is made by a three-way partition algorithm
}
\usage{
partition3wayLRMconv(lisOfDists, convMethod, single, headTrunc, tailTrunc, regridMethod = "lr", N = 256L, forDeltaCompare = 256L, useFFT = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lisOfDists}{
A list of distributions. Each distribution is a 2-column data frame
}
  \item{convMethod}{
A string indicating the convolution method.

"splitAtom4products": call convSplitAtom4products() inside

"splitAtom9products": call convSplitAtom9products() inside

"splitAtom4productsAtomUntouched": call convSplitAtom4productsAtomUntouched() inside

"splitAtom9productsAtomUntouched": call convSplitAtom9productsAtomUntouched() inside

"splitAtom4productsIrregular": call convSplitAtom4productsIrregular() inside

"splitAtom9productsIrregular": call convSplitAtom9productsIrregular() inside

"convbt": call brute-force inside
}
  \item{single}{
A boolean value indicating the truncation method. single=1 uses truncation method 1
}
  \item{headTrunc}{
Head truncation threshold
}
  \item{tailTrunc}{
Tail truncation threshold
}
  \item{regridMethod}{
A string indicating the regridding method.

"lr": linear regrid

"r3": 3 dimensional regrid

"r4": 4 dimensional regrid

"r4g": 4 dimensional regrid with Gaussian kernal
}
  \item{N}{
Maximum number of grid points allowed in the convolution distribution
}
  \item{forDeltaCompare}{
The number used to compute the delta of the convolution. Should not be less than N
}
\item{useFFT}{
Whether use FFT to convolve}
}
\details{
Please refer to the instruction of convolution for more details
}
\value{
A 2-column data frame as the convolution distribution
}
