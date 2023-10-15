\name{PMFlist2comparableHist}
\alias{PMFlist2comparableHist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert PMFs to histograms
}
\description{
Given a list of PMFs, convert it to a list of histograms that are comparable for plotting.
}
\usage{
PMFlist2comparableHist(
  distlist,
  PminmaxRatio = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distlist}{
A list of PMFs. Each PMF is a 2-column data frame. The 1st column is the value vector and the 2nd column is the probability vector.
}
  \item{PminmaxRatio}{
In each PMF, ignore the rows whose probabilities are less than the maximum probability times \code{PminmaxRatio}.
}
}


\value{
A list of histograms.
}

