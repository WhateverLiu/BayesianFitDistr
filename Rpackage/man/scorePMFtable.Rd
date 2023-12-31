% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scorePMFtable.R
\name{scorePMFtable}
\alias{scorePMFtable}
\title{Score PMF Table}
\usage{
scorePMFtable(pmfTable, claimsData)
}
\arguments{
\item{pmfTable}{PMF table as a numeric matrix. Every column is a 
(MDR, max, P0, P1, ..., Pmax).}

\item{claimsData}{Claims data as a 2-column data frame. The first column
is MDR, and the second is claim damage ratio. The MDR should be 
\eqn{\in(0, 1)}. Claim damage ratio should be \eqn{\in[0, 1]}.}
}
\value{
A list of 2 objects:
\itemize{
  \item{\code{pmfsAndScores}}{ 
    A list. The size is the number of PMFs that have been scored. 
    \code{pmfsAndScores[[i]]} is a list of 3 objects:
    \itemize{
      \item{\code{scoredPMF}}{ 
        The scored PMF as a list of 2 numeric vectors: 
        support and probability.
      }
      \item{\code{claims}}{ The claim damage ratios used to score this PMF.}
      \item{\code{ignoranceScores}}{ Ignorance score from each damage ratio.}
    }
  }
  \item{\code{nllh}}{  Negative log-likelihood of the claims data. }
}
}
\description{
Given claims data, assess the goodness of fit of the PMF table by computing
(i) ignorance scores and (ii) negative log-likelihood of the claims data.
}
\details{
For computing the negative log-likelihood of the claims data, the
PMF table is seen as a 2D joint PMF, and claims data (MDR, claim damage ratio)
are realizations of the joint. More details can be found in
\href{../doc/slides.pdf}{\code{vignette("slides", package = "NGFMfitDistr")}}
}
