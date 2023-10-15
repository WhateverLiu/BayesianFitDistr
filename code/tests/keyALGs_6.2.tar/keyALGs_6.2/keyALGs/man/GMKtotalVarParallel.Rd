\name{GMKtotalVarParallel}
\alias{GMKtotalVarParallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the covariance matrix sum given a Gaussian mixture model.
}
\description{
Compute the covariance matrix sum given a Gaussian mixture model.
}
\usage{
GMKtotalVarParallel(
  location, 
  GMkernel, 
  locationVar, 
  NofCore, 
  weightEPS = 1e-12, 
  outputCorM = F, 
  useMahalanobisDforWeight = F, 
  absorbWeightInLocalCovar = T
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{location}{
A 2-column data frame. 1st column is x-coordinate, 2nd column is y-coordinate
}
  \item{GMkernel}{
The Gaussian mixture model resulted from the EM algorithm. A list of numeric vectors. Each vector is of size 6. First element is weight, 2nd - 3rd elements is kernel mean, 4th element is var1, 5th element is covar, 6th element is var2. The first element, the weight, won't be used in the calculation. Please set it as an arbitrary value if your GM model doesn't have the mixture weights.
}
  \item{locationVar}{
A numeric vector. The variances of the location loss distributions.
}
  \item{NofCore}{
Number of threads to put in use.
}
  \item{weightEPS}{
When formulating the Gaussian mixture kernel function for a location, if the weight ends up less than weightEPS, then the Gaussian component corresponding to this weight will be deserted. weightEPS is designed to exclude the trivial components thus to control the size of the Gaussian mixture kernel function, leading to faster processing.
}
  \item{outputCorM}{
If outputCorM is true, a correlation matrix will be output.
}
\item{useMahalanobisDforWeight}{
When formulating the kernel function for a location as a mixture of the embedded kernels, should the weights be calculated depending on the probability densities the embedded kernels offer at this location or be calculated depending on the Mahalanobis distance between the location and the kernel means.
}
  \item{absorbWeightInLocalCovar}{
If absorbWeightInLocalCovar is true, a location's kernel function will be formulated as a single Gaussian function centered at this location with variance equal to that of the classic Gaussian mixture kernel function for this location.
}
}
\details{
Please refer to the non-stationary covariance analysis chapter in the Validation project report for more information
}
\value{
If outputCorM false, value of the covariance matrix sum, otherwise a list where the first element is the covariance matrix sum and the second element is the correlation matrix.
}




