\name{NGFMfitDistr}
\alias{NGFMfitDistr}
\alias{NGFMfitDistr-package}
\title{NGFMfitDistr: NGFM Distribution Fitting Toolbox}


\description{
Fit and discretize NGFM (New Generation Financial Module) compliant inflated Transformed Beta distributions. Users are strongly recommended to run the code example in \code{vignette('workflow', package = "NGFMfitDistr")} for an overview of the fitting and discretization methodologies. A more technical presentation is in \code{vignette('slides', package = "NGFMfitDistr")}.
}


\section{Fitting: from claims data to inflated TrB parameters}{
  \itemize{
    \item\code{\link{fullFit}}: An end-to-end inflated Transformed Beta fitting Pipeline.
  }
}


\section{Discretization: from inflated TrB parameters to PMF table}{
  \itemize{
    \item\code{\link{fullDiscretize}}: Full pipe for inflated Transformed Beta Discretization.
  }
}


\section{Major components in \code{\link{fullFit}}}{
  \itemize{
    \item\code{\link{computeWindows}}: Compute sliding windows.
    \item\code{\link{ensembleP0}}: Estimate \eqn{P_0}.
    \item\code{\link{estimateEmpPMFs}}: Compute empirical PMFs.
    \item\code{\link{LBFGSBtrbFitList}}: Fit TrBs.
  }
}


\section{Major components in \code{\link{fullDiscretize}}}{
  \itemize{
    \item\code{\link{FTDFQA}}: Inflated Transformed Beta Discretization.
  }
}


\section{Determine max}{
  \itemize{
    \item\code{\link{makeMax}}: Estimate proper maxes of the supports for TrB discretization.
  }
}


\section{Visualization}{
All the visualization functions are called either inside \code{\link{fullFit}} or \code{\link{fullDiscretize}}.
  \itemize{
    \item\code{\link{plotP0model}}: Plot \eqn{P_0}.
    \item\code{\link{plotTrBparameters}}: Plot TrB Parameters.
    \item\code{\link{plot45pmfsInTable}}: Plot 45 PMFs in PMF table.
    \item\code{\link{plot207pmfsInTable}}: Plot 207 PMFs in PMF table.
    \item\code{\link{plotQAparams}}: Visualize Parameter Trends in PMF table.
  }
}


\section{Miscellaneous}{
All the visualization functions are called either inside \code{\link{fullFit}} or \code{\link{fullDiscretize}}.
  \itemize{
    \item\code{\link{actuar_qtrbeta}}: Plot \eqn{P_0}.
    \item\code{\link{CharliePara}}: Naive distributed computing.
    \item\code{\link{equalChunk}}: Group jobs.
    \item\code{\link{extractMain}}: Extract main part.
    \item\code{\link{igscore}}: Ignorance Score.
    \item\code{\link{longestNonDecreasingSubseq}}: Longest Nondecreasing Subsequence.
    \item\code{\link{scorePMFtable}}: Score PMF Table.
    \item\code{\link{solve_d}}: Solve the scale parameter d.
  }
}





