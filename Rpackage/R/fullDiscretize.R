



#'
#' Full Pipe for Inflated Transformed Beta Discretization
#' 
#' Discretize inflated TrB distributions characterized by a parameter table 
#' under QA requirements. This function amalgamates \code{\link{FTDFQA}}, 
#' \code{\link{plotQAparams}}, \code{\link{plot45pmfsInTable}}, 
#' \code{\link{plot207pmfsInTable}} and \code{\link{scorePMFtable}}.
#' 
#' @inheritParams FTDFQA
#' 
#' @inheritParams fullFit
#' 
#' @param empDistrsLists  The list of empirical PMFs to which inflated TrBs 
#' are fitted. \code{empDistrsLists} is a member in the output from \code{\link{fullFit}}. 
#' The object is only needed for visualization, and would be ignored if
#' \code{figureDir=""}. Default to \code{NULL} which silences 
#' visualization.
#' 
#' @param mdrRangeInData  A numeric vector of size 2. 
#' MDR range in the original claims data. \code{mdrRangeInData}
#' is a member in the output from \code{\link{fullFit}}. The object is only needed 
#' for visualization, and would be ignored if \code{figureDir=""}. Default to
#' \code{NULL} which silences visualization.
#' 
#' @param claimsDataForScoring  Claims data for scoring the tuned PMF table.
#' See argument \code{claimsData} in \link{scorePMFtable} for more
#' details. Default to \code{NULL} which also silences scoring. 
#' 
#' @return A list of the results from \code{\link{FTDFQA}} and the results from 
#' \code{\link{scorePMFtable}}. The latter exists only if claims data are given to 
#' score the PMF table.
#'
#'
fullDiscretize = function(
    TrBtable, 
    supportSizes = c(42L, 64L), 
    regridMethod = "lr",
    RIBlib = "Numerical Recipes",
    outputProbsInRows = FALSE,
    fineDiscretizationSize = 2000,
    maxCore = 1000, 
    verbose = TRUE,
    downScaleSupport = FALSE,
    figureDir = "../figure",
    empDistrsLists = NULL,
    mdrRangeInData = NULL,
    claimsDataForScoring = NULL
    )
{
  
  rst = FTDFQA(
    TrBtable, 
    supportSizes, 
    regridMethod,
    RIBlib,
    outputProbsInRows = FALSE,
    fineDiscretizationSize,
    maxCore,
    verbose, 
    downScaleSupport)
  
  
  if (figureDir != "")
  {
    
    if (verbose) cat("Making plots...\n")
    
    
    dir.create(figureDir, showWarnings = F)
    X = rst$untuned
    
    
    Xnames = names(X)
    for (i in 1:length(X)) 
    {
      fname = paste0(figureDir, "/", Xnames[i],  
                     "-rawDiscretization-p0p1maxCv-MDR.png")
      plotQAparams(X[[i]], fname)
    }
    
    
    X = rst$tuned
    for (i in 1:length(X)) 
    {
      fname = paste0(figureDir, "/", Xnames[i],  
                     "-tunedDiscretization-p0p1maxCv-MDR.png")
      plotQAparams(X[[i]], fname)
    }
    
    
    if (!is.null(empDistrsLists)) 
    {
      
      for (i in 1:length(X))
      {
        figPath = paste0(figureDir, "/", Xnames[i],  "-pmfTable-45.pdf")
        plot45pmfsInTable(figPath, X[[i]], empDistrsLists, 
                          mdrRangeInData = mdrRangeInData)  
      }
      
      
      for (i in 1:length(X))
      {
        figPath = paste0(figureDir, "/", Xnames[i],  "-pmfTable-207.pdf")
        plot207pmfsInTable(figPath, X[[i]], empDistrsLists, 
                           mdrRangeInData = mdrRangeInData)
      }
    }
  }
  
  
  if (!is.null(claimsDataForScoring))
  {
    scores = lapply(rst$tuned, function(x) 
      scorePMFtable(x, claimsDataForScoring))
    names(scores) = names(rst$tuned)
    rst$scores = scores 
  }
  
  
  if (outputProbsInRows)
  {
    for (i in 1:length(rst))
    {
      for (j in 1:length(rst[[i]]))
        rst[[i]][[j]] = t(rst[[i]][[j]])
    }
  }

  
  rst
}







