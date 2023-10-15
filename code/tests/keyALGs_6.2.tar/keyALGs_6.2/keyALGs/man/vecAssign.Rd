\name{vecAssign}
\alias{vecAssign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
In-place modify vector
}
\description{
Alter elements of a vector without making a copy of it.
}
\usage{
vecAssign(orginalV, index, newValueV)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orginalV}{
The vector whose elements are to be replaced.
}
  \item{index}{
One-based indices of the elements to be replaced.
}
  \item{newValueV}{
Replacement values. Size should equal \code{index}'s size.
}
}
\details{
\code{orginalV[index] = newValueV}
}
\value{
None
}
