\name{makeLossDistrFromDtable}
\alias{makeLossDistrFromDtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create loss distributions
}
\description{
Given mean damage ratios, replacement values and NGFM distribution table, make a list of loss distributions.
}
\usage{
makeLossDistrFromDtable(distrTable, givenMDRs, repVals)
}
\arguments{
  \item{distrTable}{
A DATA FRAME of NGFM distribution table. A column is a vector of (MDR, max, P_1, ..., P_N) where N is the support size.
}
  \item{givenMDRs}{
A numeric vector of mean damage ratios.
}
  \item{repVals}{
A numeric vector of replacement values.
}
}
\details{
Degenerate distributions will be a data frame of 1 row.
}
\value{
A list of distributions. Each list element is a 2-column data frame.
}
