


#' Model P0
#' 
#' Model P0 as a function of MDR via the ensemble approach.
#' 
#' @param dat  A two-column data frame. The first column is the MDR. The second
#' is the claim damage ratio. Rows of the data MUST have been ordered by MDRs 
#' first and then the claim data ratios.
#'
#' @param windowSize  A positive integer. Should be less than the number of rows of
#' \code{dat}.
#' 
#' @param  slidingSpeed A positive integer. Should be less than the number of rows of
#' \code{dat}.
#' 
#' @param sampleSize  A positive integer. Sample size to create the ensemble.
#' 
#' @param NsampleSets  A positive integer. Ensemble size. Default to 100.
#' 
#' @param MDRs  A numeric vector of prescribed MDRs, usually (1e-5, 2e-5, ..., 1e-1, 
#' 1e-1 + 1e-4, 1e-1 + 2e-4, ..., 9999e-4). Default to this sequence.
#' 
#' @param randomSeed  Random seed for resampling. Default to 42.
#' 
#' @param maxCore  Maximum number of CPU cores to be used. Default to the 
#' machine's maximum number of logical processors.
#' 
#' @param tempDir  Temporary directory for saving computing progress. Default to
#' \code{'../tempFiles'} which will be created if nonexistent.
#' 
#' @param figureDir  Directory for saving plot(s). Default to
#' \code{'../figure'} which will be created if nonexistent.
#' 
#' @details The function (i) creates \code{NsampleSets} of random subsets, (ii)
#' fits \eqn{P_0} in each subset, (iii) takes the mean of \eqn{P_0}s as the 
#' final model. More details can be found in this 
#' \href{../inst/docs/slides.pdf}{slide documentation}} from Page 26.
#' 
#' @return A list of size two:
#' 
#' \code{$p0models}: A list that contains the \eqn{P_0} sequences from 
#' all the sample sets.
#' 
#' \code{$P0}: Mean of all the \eqn{P_0} sequences in \code{$p0models}.
#' 
#' @example inst/examples/ensembleP0.R
ensembleP0 = function(
  dat,
  windowSize,
  slidingSpeed,
  sampleSize,
  NsampleSets = 100L,
  MDRs = c(1:10000, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5,
  randomSeed = 42L,
  maxCore = parallel::detectCores(),
  tempDir = "../tempFiles/CharlieTempMP/C",
  figureDir = "../figure"
  )
{
  targetMDRs = MDRs
  colnames(dat) = c("MDR", "CDR")
  commonData = list(datResv = dat, windowSize = windowSize, 
                    slidingSpeed = slidingSpeed, sampleSize = sampleSize, 
                    targetMDRs = MDRs)
  f = function(seed, commonData)
  {
    datResv = commonData$datResv
    windowSize = commonData$windowSize
    slidingSpeed = commonData$slidingSpeed
    sampleSize = commonData$sampleSize
    targetMDRs = commonData$targetMDRs
    dat = datResv
    set.seed(seed)
    dat = dat[sample(nrow(dat), sampleSize), , drop = FALSE]
    dat = dat[order(dat$MDR), , drop = FALSE]
    windows = computeWindows(sampleSize, windowSize, 
                             start = 1L, speed = slidingSpeed)
    modelP0(
      mdr = dat$MDR, cdr = dat$CDR, windows = windows, 
      targetMDRs = targetMDRs,
      eps = 1e-10, P0modeltype = "hyman", warnUsers = FALSE)
  }
  
  
  cat("Computing the ensemble of P0...\n\n")
  
  
  tempDir = paste0(tempDir, "/CharlieTempMP/C")
  p0models = CharliePara(
    X = as.list(seq(rseed, len = NsampleSets, by = 1L)),
    commonData = commonData, fun = f, maxNprocess = maxCore, 
    MPtempDir = tempDir)
  
  
  cat("\n")
  
  
  P0 = rowMeans(as.data.frame(p0models))
  
  
  if (!is.null(figureDir))
  {
    dir.create(figureDir, showWarnings = FALSE)
    cat("Making plot...\n")  
    fname = paste0(figureDir, "/p0modelsFigure.png")
    png(fname, width = 10, height = 0.618 * 10, 
        unit = "in", res = 120)
    par(mar = c(4, 5, 0, 0), family = "serif")
    for (i in 1:length(p0models))
    {
      if (i == 1) plt = plot else plt = lines
      plt(x = targetMDRs, y = p0models[[i]], type = "l", col = "gray",
          xlab = "MDR", ylab = expression(P[0]), cex.lab = 2, 
          cex.axis = 1.5, las = 1, bty = "L")
    }
    lines(x = targetMDRs, y = P0, col = "black", lwd = 2)
    legend("topright", legend = c(
      "Hyman interpolation\nensemble, size = 100\n", "Mean"), 
      bty = "n", cex = 2, lwd = c(2, 2), col = c("gray", "black"))
    dev.off()
  }
  
  
  cat("Plot has been saved as ", fname, "\n\n")
  
  
  list(p0models = p0models, P0 = P0)
}