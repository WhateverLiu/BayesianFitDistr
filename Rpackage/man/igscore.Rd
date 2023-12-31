% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{igscore}
\alias{igscore}
\title{Ignorance Score}
\usage{
igscore(X, x, redistributePwithin = TRUE)
}
\arguments{
\item{X}{A PMF as a list or dataframe of two numeric vectors. 
The first is the support. The second is the probabilities.}

\item{x}{A vector of realizations.}

\item{redistributePwithin}{A boolean. When a realization falls within the
range of the PMF's support, \code{redistributePwithin = FALSE} will select
the nearest probability mass for calculation, \code{TRUE} will interpolate
the two nearest probability masses for calculation. Default to \code{TRUE}.}
}
\value{
A numeric vector of the ignorance scores.
}
\description{
Compute ignorance score given a PMF and a vector of realizations.
}
\details{
For a single realization \eqn{x}, the ignorance score = 
\eqn{-\log_2 f(x)} where \eqn{f} is the PDF. More details can be found
\href{https://journals.ametsoc.org/view/journals/mwre/130/6/1520-0493_2002_130_1653_epfuit_2.0.co_2.xml?tab_body=fulltext-display}{here}.
}
