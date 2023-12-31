% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CharlieMPlib.R
\name{equalChunk}
\alias{equalChunk}
\title{Group jobs}
\usage{
equalChunk(ns, Nchunks)
}
\arguments{
\item{ns}{A numeric or integer vector. \code{ns[i]} is the computing cost
for the i-th job.}

\item{Nchunks}{The number of groups for the jobs.}
}
\value{
A list of index vectors.
}
\description{
Group jobs such that the computing cost for each group is approximately the 
same.
}
