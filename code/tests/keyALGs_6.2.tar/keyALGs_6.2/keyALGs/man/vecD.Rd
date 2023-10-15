\name{vecD}
\alias{vecD}
\title{vecD}
\description{
Compute the distance between two vectors
}
\usage{
vecD(x, y, ms = "Cheb")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{A numeric vector.}
\item{y}{A numeric vector of the same size of \code{x}}
\item{ms}{A character string. Name of the distance.}
}
\details{The function will do a size check. If the x and y have different length, the function will return an error message.

Please refer to \url{http://citeseerx.ist.psu.edu/viewdoc/download?rep=rep1&type=pdf&doi=10.1.1.154.8446} to see definitions of the distances below.

\code{ms = "Cheb"} returns the Komogorov distance. Definition 4

\code{ms = "Euc"} returns the Euclidean distance. Definition 1

\code{ms = "CB"} returns the City Block distance. Definition 2

\code{ms = "Can"} returns the Canberra distance. Definition 9

\code{ms = "B"} returns the Bhattacharyya distance. Definition 25

\code{ms = "H"} returns the Hellinger distance. Definition 26

\code{ms = "J"} returns the Jeffreys distance. Definition 38

\code{ms = "Top"} returns the Topsoe distance. Definition 40

\code{ms = "Tani"} returns the  Tanimoto distance. Definition 17
}

\value{A numeric value as the distance.}


