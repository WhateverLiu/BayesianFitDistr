\name{keyALGs}
\alias{keyALGs}
\alias{keyALGs-package}
\title{keyALGs: Key Algorithms}


\description{
Key algorithms for financial modeling research. Many are for exploration and not introduced to the software team. Some are helper functions for the author's personal use. New users are recommended to start with functions listed on this page. They are the most frequently used, and mainly deal with PMFs (probability mass functions).
}


\section{Regriding}{
  \itemize{
    \item\code{\link{rglr}}: Linear regriding. Closest to the regriding routine in software release.
    \item\code{\link{rglr4split}}: Four-point regriding.
    \item\code{\link{rglr4splitDetails}}: Four-point regriding with process information returned.
    \item\code{\link{LMMdiscreteToDiscrete}}: Local moment matching.
  }
}


\section{Empirical PMF from samples}{
  \itemize{
    \item\code{\link{emp}}: Generate PMF that characterizes distribution of a sample set via regriding. Function name stands for "empirical".
  }
}


\section{Discretization}{
  \itemize{
    \item\code{\link{LMMcontinousToDiscrete}}: Local moment matching for discretization of a continous distribution.
  }
}


\section{Convolve two distributions}{
  \itemize{
    \item\code{\link{convSplitAtom4productsIrregular}}: Four-product convolution. Closest to the convolution routine in software release. WARNING!!: Not to be confused with \code{\link{convSplitAtom4products}}, which is an experimental routine that has not been used since 2016.
    \item\code{\link{convBt}}: Brute-force convolution. Able to convolve distributions with fully irregular supports.
  }
}


\section{Convolve many distributions}{
  \itemize{
    \item\code{\link{sequentialConvAll}}: Sequentially convolve a list of distributions.
    \item\code{\link{sortedSequentialConvAll}}: Sort the list of distributions by distribution's support step size, and then call \code{\link{sequentialConvAll}()}. This MAY mitigate variance inflation.
    \item\code{\link{closestPairConvAll}}: Recursively searches for pairs of PMFs that have the closest step sizes to convolve. This typically incurs the least variance inflation. It is later also referred to as \emph{recursive nearest neighbor order of convolution}.
  }
In the software release, the order of convolution also depends on the risk groupings by financial terms or by geographic regions, etc.
}


\section{Comonotonize two distributions}{
  \itemize{
    \item\code{\link{como}}: Compute the sum of two distributions assuming maximum correlation.
  }
}


\section{Comonotonize many distributions}{
  \itemize{
    \item\code{\link{comoDistList}}: Comonotonize a list of distributions in sequential order.
  }
}


\section{Counter-comonotonize two distributions}{
  \itemize{
    \item\code{\link{counterComo}}: Compute the sum of two distributions assuming minimum correlation.
  }
}


\section{Distribution mean and variance}{
  \itemize{
    \item\code{\link{Mean}}: Compute the mean of a PMF.
    \item\code{\link{Var}}: Compute the variance of a PMF.
  }
}
  
  
\section{Max-convolution, Aggregate PMFs with mixture weights, correlations or copulas}{
  \itemize{
    \item See package \code{\link{keyALGsExt}}. NGFM research algorithms are still in development and have adopted newer C++ standards. These algorithms deal with a quite different set of problems, thus a separate package.
  }
}

















