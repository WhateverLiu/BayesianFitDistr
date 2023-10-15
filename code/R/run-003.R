setwd('/finance_develop/Charlie/BayesianFitDistr/code')


source("R/rfuns.R")
source("R/CharlieMPlib.R")
Rcpp::sourceCpp("src/amalgamation.cpp", verbose = T, cacheDir = "../tempFiles")


# ==============================================================================
# Read data in.
# ==============================================================================
oldDistTable = data.table::setDF(data.table::fread(
  "../tempFiles/CovA_nPts64-debora.csv")) # For comparing ignorance score.
load("../data/cvgAdata.Rdata") # Coverage A
# load("../data/cvgCdata.Rdata") # Coverage C


# A claim’s damage ratio is considered to be zero if it is <= 0.002


oldDistTable = t(oldDistTable)
dimnames(oldDistTable)= NULL


dat = cvgData
dat$CDR[dat$CDR <= 0.002] = 0
datResv = dat[order(dat$MDR, dat$CDR), , drop = F]
rm(cvgData); gc()
# ==============================================================================




# ==============================================================================
# Read in old distribution table for Bayesian update
# ==============================================================================
if (F)
{
  tohoku = data.table::setDF(data.table::fread("tests/ResCovA64_Final.csv"))[-1]
  tohoku = t(tohoku)
  dimnames(tohoku) = NULL
  tohoku = apply(tohoku, 2, function(x)
  {
    list(val = seq(0, x[2], len = length(x) - 2L), P = x[-c(1, 2)])
  })
  tohoku = tohoku[-c(1, length(tohoku))]  
}




# ==============================================================================
# Sliding window size and speed.
# ==============================================================================
windowSizeRatio = 0.01
slidingSpeedRatioAsWindowSize = 0.01
maxQuantile = 0.99
maxLowerBound = 0.05




# ==============================================================================
# Model P0s, direct optimization for ignorance score.
# ==============================================================================
if (F)
{
  
  dat = datResv[c("MDR", "CDR")]
  tmpGrouped = aggregate(dat["CDR"], list(
    mdrInd = as.integer(round(dat$MDR * 1e5))), function(x) x)
  tmp = tmpGrouped
  tmp = data.frame(
    MDR = tmp$mdrInd / 1e5, 
    p0 = unlist(lapply(
      tmp$CDR, function(x) sum(x < 1e-10) / length(x))), 
    Nzero = 
      unlist(lapply(tmp$CDR, function(x) sum(x < 1e-10))) )
  tmp$logP0 = log(tmp$p0)
  dat = tmp
  
  
  # ============================================================================
  # Assign colors to data points to be ploted.
  # ============================================================================
  nzeroRanks = rank(dat$Nzero, ties.method = "average")
  tmp = qnorm((1:length(nzeroRanks)) / (length(nzeroRanks) + 1))
  tmp = tmp[nzeroRanks]
  color = colorRampPalette(c("white", "gray", "black"))(length(tmp))
  whichColor = as.integer(round((tmp - min(tmp)) / diff(range(tmp)) * (
    length(tmp) - 1) + 1))
  
  
  # ============================================================================
  # Plot MDR vs logP0 with color
  # ============================================================================
  tmpDat = dat[c("MDR", "logP0", "Nzero")]
  tmpDat = tmpDat[is.finite(tmpDat$logP0), , drop = F]
  cutPercentage = 0.99
  mdrCut = datResv$MDR[round(nrow(datResv) * cutPercentage)]
  datForMonopoly = tmpDat[tmpDat$MDR <= mdrCut, , drop = F]
  monopoly = MonoPoly::monpol(logP0 ~ MDR, data = datForMonopoly, degree = 5,
                              a = datForMonopoly$MDR[1], b = mdrCut, 
                              weights = datForMonopoly$Nzero )
  fittedLogP0 = fitted(monopoly)
  names(fittedLogP0) = NULL
  targetMDRs = c(1:1e4L, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5
  ind = which(targetMDRs >= tmpDat$MDR[1] & targetMDRs <= mdrCut)
  inRangePred = predict(monopoly, newdata = list(MDR = targetMDRs[ind]))
  attributes(inRangePred) = NULL
  
  
  headMDRs = c(0, targetMDRs[ind[1:2]])
  headP0s = c(1, exp(inRangePred[1:2]))
  cbfun = splinefun(x = headMDRs, y = headP0s, method = 'hyman')
  headP0s = cbfun(targetMDRs[1:(ind[1] - 1L)])
  
  
  tailMDRs = c(targetMDRs[ind[c(length(ind)-1, length(ind))]], 1)
  tailP0s = c(exp(inRangePred[c(length(ind)-1, length(ind))]), 0)
  cbfun = splinefun(x = tailMDRs, y = tailP0s, method = 'hyman')
  tailP0s = cbfun(targetMDRs[(ind[length(ind)] + 1L):length(targetMDRs)])
  
  
  P0 = c(headP0s, exp(inRangePred), tailP0s)
  sum(log2(P0[match(tmpGrouped$mdrInd, as.integer(round(targetMDRs * 1e5)))]) * 
    unlist(lapply(tmpGrouped$CDR, function(x) sum(x < 1e-10)))) / nrow(datResv)
  
  
  
  
  png("../figure/MDRvsLogP0.png", width = 12, height = 6, unit = "in", 
      res = 240)
  par(mar = c(4.2, 5, 1, 1), family = "serif", mfrow = c(1, 2))
  plot(tmpDat$MDR, tmpDat$logP0, cex = 0.5, pch = 16, col = color[whichColor], 
       bty = "L", xlab = "MDR", ylab = expression(logP[0]), cex.lab = 2,
       cex.axis = 1.5, las = 1)
  lines(datForMonopoly$MDR, fittedLogP0, col = "blue")
  lines(x = c(mdrCut, mdrCut), y = c(-1e100, 1e100), type = "h", lty = 2)
  legend("right", legend = paste0(
    "Fit the\nlowest ", round(cutPercentage * 100, 2), 
    "%\n MDRs only" ), bty = "n", cex = 1.5)
  
  
  plot(tmpDat$MDR, exp(tmpDat$logP0), cex = 0.5, pch = 16, 
       col = color[whichColor], 
       bty = "L", xlab = "", ylab = expression(P[0]), cex.lab = 2,
       cex.axis = 1.5, las = 1, xlim = c(0, 1))
  lines(targetMDRs, P0, col = "blue")
  lines(c(0, targetMDRs[1]), c(1, P0[1]), col = "blue")
  lines(c(targetMDRs[length(targetMDRs)], 1), 
        c(P0[length(P0)], 0), col = "blue")
  dev.off()
  
}




# ==============================================================================
# Model P0s, ensemble
# ==============================================================================
if (T)
{
  
  
  sampleSize = as.integer(round(nrow(datResv) * (1-1/exp(1)) ))
  NsampleSets = 100L
  rseed = 42L
  
  
  windowSize = max(1L, as.integer(round(sampleSize * windowSizeRatio)))
  slidingSpeed = as.integer(max(1, round(
    windowSize * slidingSpeedRatioAsWindowSize)))
  windows = computeWindows(sampleSize, windowSize, 
                           start = 1L, speed = slidingSpeed)
  tmp = datResv$MDR
  mdrRangeInData = c(mean(tmp[1:(windowSize - 1L)]), 
                     mean(tmp[(length(tmp) - windowSize + 1L):length(tmp)]))
  
  
  targetMDRs = c(1:1e4L, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5
  
  
  # ============================================================================
  # Parallel run
  # ============================================================================
  if (T)
  {
    
    
    commonData = list(datResv = datResv, windowSize = windowSize, 
                      slidingSpeed = slidingSpeed, sampleSize = sampleSize, 
                      targetMDRs = targetMDRs)
    f = function(seed, commonData)
    {
      datResv = commonData$datResv
      windowSize = commonData$windowSize
      slidingSpeed = commonData$slidingSpeed
      sampleSize = commonData$sampleSize
      targetMDRs = commonData$targetMDRs
      dat = datResv
      set.seed(seed)
      dat = dat[sample(nrow(dat), sampleSize), , drop = F]
      dat = dat[order(dat$MDR), , drop = F]
      windows = computeWindows(sampleSize, windowSize, 
                               start = 1L, speed = slidingSpeed)
      modelP0(
        mdr = dat$MDR, cdr = dat$CDR, windows = windows, 
        targetMDRs = targetMDRs,
        eps = 1e-10, P0modeltype = "hyman", warnUsers = FALSE)
    }
    
    
    p0models = CharliePara(X = as.list(seq(rseed, len = NsampleSets, by = 1L)), 
                           commonData = commonData, fun = f, maxNprocess = 15L)
    
    
  } else # Single core run.
  {
    
    
    p0models = list()
    for (i in 1:NsampleSets)
    {
      cat(i, "")
      set.seed(i + rseed - 1L) # To be consistent with the Parallel solution.
      dat = datResv[c("MDR", "CDR")]
      dat = dat[sample(nrow(dat), sampleSize, replace = F), , drop = F]
      dat = dat[order(dat$MDR), , drop = F]
      p0models[[length(p0models) + 1]] = modelP0(
        mdr = dat$MDR, cdr = dat$CDR, windows = windows, 
        targetMDRs = targetMDRs,
        eps = 1e-10, P0modeltype = "hyman", warnUsers = FALSE)
    }
  }
  
  
  P0 = rowMeans(as.data.frame(p0models))
  
  
  # Plot the models.
  if (F)
  {
    
    png("../figure/p0modelsFigure.png", width = 10, height = 0.618 * 10, 
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
  
  
  if (!all(P0 < 1 - targetMDRs))
    stop("P0s are too high: P0 >= 1 - targetMDR. ",
         "Only Bernoulli distributions could be fitted to the data. ",
         "Remove some zero claims to make the Transformed Beta model feasible.")
  if (F)
  {
    plot(x = targetMDRs, y = P0, type = "l")  
  }
  
  
  
}




# ==============================================================================
# Model P0s, deterministic
# ==============================================================================
if (F)
{
  
  
  # sampleSize = as.integer(round(nrow(datResv) * (2/3)))
  sampleSize = nrow(datResv)
  NsampleSets = 1L
  windowSize = max(1L, as.integer(round(sampleSize * windowSizeRatio)))
  slidingSpeed = as.integer(max(1, round(
    windowSize * slidingSpeedRatioAsWindowSize)))
  windows = computeWindows(sampleSize, windowSize, 
                           start = 1L, speed = slidingSpeed)
  tmp = datResv$MDR
  mdrRangeInData = c(mean(tmp[1:(windowSize - 1L)]), 
                     mean(tmp[(length(tmp) - windowSize + 1L):length(tmp)]))
  
  
  targetMDRs = c(1:1e4L, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5
  p0models = list()
  set.seed(123)
  for (i in 1:NsampleSets)
  {
    cat(i, "")
    dat = datResv[c("MDR", "CDR")]
    dat = dat[sample(nrow(dat), sampleSize, replace = F), , drop = F]
    dat = dat[order(dat$MDR), , drop = F]
    p0models[[length(p0models) + 1]] = modelP0(
      mdr = dat$MDR, cdr = dat$CDR, windows = windows, 
      targetMDRs = targetMDRs,
      eps = 1e-10, P0modeltype = "hyman", warnUsers = FALSE)
  }
  P0 = rowMeans(as.data.frame(p0models))
  if (!all(P0 < 1 - targetMDRs))
    stop("P0s are too high: P0 >= 1 - targetMDR. ",
         "Only Bernoulli distributions could be fitted to the data. ",
         "Remove some zero claims to make the Transformed Beta model feasible.")
  plot(x = targetMDRs, y = P0, type = "l")
  
  
}




# ==============================================================================
# Model P0s, round and interpolate
# ==============================================================================
if (F) # Not working. Fluctuates crazy.
{
  
  dat = datResv[c("MDR", "CDR")]
  tmp = aggregate(dat["CDR"], list(
    mdrInd = as.integer(round(dat$MDR * 1e5))), function(x) x)
  tmp = data.frame(MDR = tmp$mdrInd / 1e5, p0 = unlist(lapply(
    tmp$CDR, function(x) sum(x < 1e-10) / length(x))))
  
  
}




# ==============================================================================
# Make emp PMFs, run bidirectional sequential fitting. No ensemble.
# ==============================================================================
if (T)
{
  
  source("R/rfuns.R")
  datNonzeroResv = datResv[datResv$CDR > 1e-10, c("MDR", "CDR")]
  datNonzeroResv = datNonzeroResv[order(
    datNonzeroResv$MDR, datNonzeroResv$CDR), , drop = F]
  windowSize = max(1L, as.integer(round(nrow(datNonzeroResv) * windowSizeRatio)))
  slidingSpeed = as.integer(max(1, round(
    windowSize * slidingSpeedRatioAsWindowSize)))
  
  
  abcLB = c(1.01, 0.1, 0.1)
  abcUB = c(30, 30, 30)
  abc = matrix(c(4, 5, 6))
  
  
  dat = datNonzeroResv
  set.seed(456)
  dat = dat[sample(nrow(dat)), , drop = F]
  dat = dat[order(dat$MDR), , drop = F]
  
  
  windows = computeWindows(nrow(dat), windowSize, start = 1L, 
                           speed = slidingSpeed)
  
  
  tmp = datResv$MDR
  mdrRangeInData = c(mean(tmp[1:(windowSize - 1L)]), 
                     mean(tmp[(length(tmp) - windowSize + 1L):length(tmp)]))
  
  
  empDistrsLists = estimateEmpPMFs(
    mdr = dat$MDR, cdr = dat$CDR, windows, P0 = P0,
    targetMDRs = targetMDRs,
    regridMethod = "lmm",
    nonDegenerateOldDistrs = NULL,
    sampleWeightOnOldDistrs = 0.5,
    empDistrSupportSize = 64L,
    maxCore = 1000L)
  
  
  # ============================================================================
  # TRUE Use bias-corrected version as the target or not.
  # ============================================================================
  if (F) empDistrs = empDistrsLists$biasCorrected else empDistrs =
    empDistrsLists$biased
  
  
  tmp = extractMain(
    empDistrs, targetMDRs, normalizeMainPart = TRUE, removeZeroPs = FALSE)
  condEmpDistrs = tmp$conditionalDistrs
  lm1 = tmp$lm1
  
  
  startMDRtoFit = median(datResv$MDR)
  startingPmf = which.min(abs(targetMDRs - startMDRtoFit))
  cat("Which PMF will be fitted first =", startingPmf, "\n")
  
  
  optRst = LBFGSBtrbFitList(
    abc = abc, 
    lm1 = lm1,
    empDistrList = condEmpDistrs,
    abcLB = abcLB, abcUB = abcUB, 
    scaleEps = 1e-8, scaleMaxit = 100, 
    distanceFun = "likelihood",
    max_iterations = 100, maxCore = 1000, 
    RIBlib = "Numerical Recipes", 
    sequentialUpdate = startingPmf, 
    hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
    epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
    max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
    ftol = 1e-4, wolfe = 0.9)
  
  
  # ============================================================================
  # Plot transition of parameters.
  # ============================================================================
  if (F)
  {
    
    # png("../figure/TrBtransitionParam-5-1-2.png", width = 9, height = 9 * 0.618,
    #     unit = "in", res = 120)
    # png("../figure/TrBtransitionParam-5-5-2.png", width = 9, height = 9 * 0.618,
    #     unit = "in", res = 120)
    # png("../figure/TrBtransitionParam-2-3-4.png", width = 9, height = 9 * 0.618,
    #     unit = "in", res = 120)
    png("../figure/TrBtransitionParam.png", width = 9, height = 9 * 0.618,
        unit = "in", res = 120)
    par(mar = c(4.2, 5, 0, 0), family ="serif")
    layout(rbind(c(1, 2),
                 c(3, 4),
                 c(5, 5)))
    plot(optRst$param[1,], xlab = "", ylab = "a", 
         cex.lab = 2, cex.axis = 1.5, type = "l", bty = "L", las = 1, 
         col = "darkblue", xaxt = "n")
    axis(side = 1, at = c(0, 5000, 10000, 15000, 18999), labels = 
           c(0, 0.05, 0.1, 0.6, 1), cex.axis = 1.5)
    lines(x = c(startingPmf, startingPmf), y = c(-1e10, 1e10), lty = 2)
    plot(optRst$param[2,], xlab = "", ylab = "b", 
         cex.lab = 2, cex.axis = 1.5, type = "l", bty = "L", las = 1, 
         col = "darkblue", xaxt = "n")
    axis(side = 1, at = c(0, 5000, 10000, 15000, 18999), labels = 
           c(0, 0.05, 0.1, 0.6, 1), cex.axis = 1.5)
    lines(x = c(startingPmf, startingPmf), y = c(-1e10, 1e10), lty = 2)
    plot(optRst$param[3,], xlab = "", ylab = "c", 
         cex.lab = 2, cex.axis = 1.5, type = "l", bty = "L", las = 1, 
         col = "darkblue", xaxt = "n")
    axis(side = 1, at = c(0, 5000, 10000, 15000, 18999), labels = 
           c(0, 0.05, 0.1, 0.6, 1), cex.axis = 1.5)
    lines(x = c(startingPmf, startingPmf), y = c(-1e10, 1e10), lty = 2)
    param = optRst$param
    d = solve_d(param, eps= 1e-8, maxit = 100)
    plot(d, xlab = "", ylab = "d", 
         cex.lab = 2, cex.axis = 1.5, type = "l", bty = "L", las = 1, 
         col = "darkblue", xaxt = "n")
    axis(side = 1, at = c(0, 5000, 10000, 15000, 18999), labels = 
           c(0, 0.05, 0.1, 0.6, 1), cex.axis = 1.5)
    lines(x = c(startingPmf, startingPmf), y = c(-1e10, 1e10), lty = 2)
    plot(optRst$fval, xlab = "MDR", ylab = "", 
         cex.lab = 2, cex.axis = 1.5, type = "l", bty = "L", las = 1, 
         col = "darkblue", xaxt = "n")
    axis(side = 1, at = c(0, 5000, 10000, 15000, 18999), labels = 
           c(0, 0.05, 0.1, 0.6, 1), cex.axis = 1.5)
    title(ylab = "Objective function", cex.lab = 1.5, line = 3.6)
    lines(x = c(startingPmf, startingPmf), y = c(-1e10, 1e10), lty = 2)
    legend("bottomright",bty = "n", legend = paste0("Starting MDR"), 
           lwd = 1, lty = 2, cex = 2)
    ycoor = max(optRst$fval) - diff(range(optRst$fval)) * 0.1
    xcoor = startingPmf
    text(x = xcoor, y = ycoor, labels = paste0("Initial a b c = ", 
         paste0(abc, collapse = " ")), cex = 2)
    dev.off()
    
    
  }
  
  
  # ============================================================================
  # Merge pngs into a pdf.
  # ============================================================================
  if (F)
  {
    
    
    magick::image_write(magick::image_read( c(
      "../figure/TrBtransitionParam.png", 
      "../figure/TrBtransitionParam-5-1-2.png", 
      "../figure/TrBtransitionParam-5-5-2.png", 
      "../figure/TrBtransitionParam-2-3-4.png") ), 
      "../figure/TrBtransitionParamMerged.pdf", 
      format = "pdf")
    
    
  }
  
  
}




# ==============================================================================
# Discretization. Round 1.
# ==============================================================================
if (T)
{
  
  param = optRst$param
  d = solve_d(param, eps= 1e-8, maxit = 100)
  abcd = rbind(param[1:3, , drop = F], d)
  tailProbThreshold = 1 - maxQuantile
  
  
  # ============================================================================
  # Make and lower bound the maxes.
  # ============================================================================
  if (F)
  {
    maxes = pmin(1, actuar_qtrbeta(
      1 - tailProbThreshold / (1 - P0), abcd))
    maxes[!is.finite(maxes)] = 0
    maxes = pmax(maxLowerBound, maxes)
    ind = longestNonDecreasingSubseq(maxes)
    maxFun = splinefun(x = targetMDRs[ind], y = maxes[ind], method = "hyman")
    maxes = maxFun(targetMDRs)
  } else # Use the old distribution's maxes.
  {
    maxes = oldDistTable[2, -c(1, ncol(oldDistTable))]
  }
  TrBtable = rbind(targetMDRs, maxes, P0, abcd)
  
  
  
  
  # ============================================================================
  # Discretize and tune PMF table.
  # ============================================================================
  pmftablelist = FTDFQA(
    TrBtable, supportSize = c(64), 
    regridMethod = "lr",
    RIBlib = "Numerical Recipes",
    outputProbsInRows = F,
    fineDiscretizationSize = 2000,
    maxCore = 100, verbose = T, downScaleSupport = F)
  pmftable = pmftablelist$distTable
  pmftableBT = pmftablelist$distTableBeforeTune
  
  
  # ============================================================================
  # Visualization.
  # ============================================================================
  if (F)
  {
    # ==========================================================================
    # Plot MDR vs P0, P1, maxes, cv for raw discretization.
    # ==========================================================================
    if (T)
    {
      
      # pdf("../figure/rawDiscretization-p0p1maxCv-MDR.pdf", width = 10, 
      #     height = 10 * 0.618)
      png("../figure/rawDiscretization-p0p1maxCv-MDR.png", width = 10, 
          height = 10 * 0.618, res = 120, unit = "in")
      par(mar = c(4.2, 5, 0, 0), family = "serif")
      layout(rbind(c(1, 1, 2, 2, 3, 3),
                   c(4, 5, 5, 6, 6, 7)))
      exclu = -c(1L, ncol(pmftableBT))
      # xs = 1:length(targetMDRs)
      xs = targetMDRs
      plot(x = xs, y = pmftableBT[2, exclu], type = "l", 
           col = "darkblue", bty = "L", xlim = range(xs), ylim = c(0, 1),
           ylab = "Max", cex.lab = 2, cex.axis = 1.5, las = 1, xlab = "" )
      plot(x = xs, y = pmftableBT[3, exclu], col = "darkblue", type = "l", 
           ylab = expression(P[0]), cex.lab = 2, cex.axis = 1.5, bty = "L", 
           las = 1, xlab = "" )
      plot(x = xs, y = pmftableBT[nrow(pmftableBT), exclu], 
           col = "darkblue", cex.lab = 2, cex.axis = 1.5, bty = "L", 
           las = 1, ylab = expression(P[max]), type = "l", xlab = "" )
      plot.new()
      ms = apply(pmftableBT[, exclu], 2, function(x)
      {
        sum(seq(0, x[2], len = nrow(pmftableBT) - 2) * x[-c(1, 2)])
      })
      plot(x = xs, y = ms, col = "darkblue", cex.lab = 2, cex.axis = 1.5, 
           bty = "L", las = 1, ylab = expression(mu), type = "l", xlab = "MDR" )
      sds = sqrt(apply(pmftableBT[, exclu], 2, function(x)
      {
        sum(seq(0, x[2], len = nrow(pmftableBT) - 2) ^ 2 * x[-c(1, 2)])
      }) - ms ^ 2)
      cv = sds / ms
      # cvScaled = (cv - min(cv)) / diff(range(cv))
      plot(x = xs, y = cv, col = "darkblue", cex.lab = 2, cex.axis = 1.5, 
           bty = "L", las = 1, ylab = expression(sigma/mu), type = "l", xlab = "" )
      # axis(side = 4, at = c(0, 0.5, 1), 
      #      labels = round(c(0, diff(range(cv)) / 2, max(cv)), 3), las = 1)
      dev.off()
    }
    
    
    # ==========================================================================
    # Plot MDR vs P0, P1, maxes, cv for tuned discretization.
    # ==========================================================================
    if (T)
    {
      png("../figure/tunedDiscretization-p0p1maxCv-MDR.png", width = 10, 
          height = 10 * 0.618, unit = "in", res = 120)
      par(mar = c(4.2, 5, 0, 0), family = "serif")
      layout(rbind(c(1, 1, 2, 2, 3, 3),
                   c(4, 5, 5, 6, 6, 7)))
      exclu = -c(1L, ncol(pmftable))
      xs = targetMDRs
      plot(x = xs, y = pmftable[2, exclu], type = "l", 
           col = "darkblue", bty = "L", xlim = range(xs), ylim = c(0, 1),
           ylab = "Max", cex.lab = 2, cex.axis = 1.5, las = 1, xlab = "" )
      plot(x = xs, y = pmftable[3, exclu], col = "darkblue", type = "l", 
           ylab = expression(P[0]), cex.lab = 2, cex.axis = 1.5, bty = "L", 
           las = 1, xlab = "" )
      plot(x = xs, y = pmftable[nrow(pmftable), exclu], 
           col = "darkblue", cex.lab = 2, cex.axis = 1.5, bty = "L", 
           las = 1, ylab = expression(P[max]), type = "l", xlab = "" )
      plot.new()
      ms = apply(pmftable[, exclu], 2, function(x)
      {
        sum(seq(0, x[2], len = nrow(pmftable) - 2) * x[-c(1, 2)])
      })
      plot(x = xs, y = ms, col = "darkblue", cex.lab = 2, cex.axis = 1.5, 
           bty = "L", las = 1, ylab = expression(mu), type = "l", xlab = "MDR" )
      sds = sqrt(apply(pmftable[, exclu], 2, function(x)
      {
        sum(seq(0, x[2], len = nrow(pmftable) - 2) ^ 2 * x[-c(1, 2)])
      }) - ms ^ 2)
      cv = sds / ms
      plot(x = xs, y = cv, col = "darkblue", cex.lab = 2, cex.axis = 1.5, 
           bty = "L", las = 1, ylab = expression(sigma/mu), type = "l", xlab = "" )
      dev.off()  
    }
    
    
    # ==========================================================================
    # Make 45 plots.
    # ==========================================================================
    if (T)
    {
      tmp = unique(as.integer(round(c(seq(1e-5, 1e-4, by = 1e-5), 
                                      seq(1e-4, 1e-3, by = 1e-4), 
                                      seq(1e-3, 1e-2, by = 1e-3), 
                                      seq(1e-2, 1e-1, by = 1e-2), 
                                      seq(1e-1, 1-1e-1, by = 1e-1)) * 1e5)))
      whichEmpDistrs = match(tmp, as.integer(round(targetMDRs * 1e5)))
      empDistrsToPlot = empDistrsLists$biased[whichEmpDistrs]
      empDistrsBiasCorrectedToPlot = empDistrsLists$biasCorrected[whichEmpDistrs]
      whichPMFinTable2plot = match(tmp, as.integer(round(pmftable[1,] * 1e5)))
      tablePmf2plot = pmftableToDistrList(pmftable[, whichPMFinTable2plot, drop = F])
      pdf("../figure/pmfTable-45.pdf", width = 10, height = 10 * 0.618)
      par(mar = c(4.2, 5, 0, 0), mfrow = c(2, 1), family = "serif")
      options(scipen = 999)
      for (i in 1:length(tablePmf2plot))
      {
        X = empDistrsToPlot[[i]]
        X = keyALGs::rglr(X, ngd = seq(
          X$val[1], X$val[length(X$val)], len = length(X$val)))
        Y = tablePmf2plot[[i]]
        Z = empDistrsBiasCorrectedToPlot[[i]]
        Z = keyALGs::rglr(Z, ngd = seq(
          Z$val[1], Z$val[length(Z$val)], len = length(Z$val)))
        Xhist = wrapAsHist(X)
        s = (diff(range(Y$val)) / (length(Y$val) - 1)  ) / 
          (diff(range(X$val)) / (length(X$val) - 1)  )
        l = length(Xhist$mids)
        if (Xhist$mids[l] >= 1-1e-10) 
          Xhist$density[-c(1, l)] = Xhist$density[-c(1, l)] * s
        else Xhist$density[-c(1)] = Xhist$density[-c(1)] * s
        xlim = c(0, max(c(range(Y$val), range(Xhist$breaks), range(Z$val))))
        ylim = range(c(0, max(Y$P), max(Xhist$density), max(Z$P)))
        
        
        hcol = "gray"
        plot(Xhist, freq = F, border = NA, col = hcol, main = "", 
             xlab = "", cex.lab = 2, cex.axis = 1.5, ylab = "P",
             xlim = xlim, ylim = ylim, las = 1)
        lines(as.data.frame(Y), type = "h", col = "black", lwd = 1)
        tmpcoor = keyALGs::plotCoor(0.5, 0.92)
        fullmdr = keyALGs::Mean(Y)
        Nround = 5
        
        
        mdr = keyALGs::num2str(keyALGs::Mean(Y), Nround = Nround)
        text(x = tmpcoor[1], y = tmpcoor[2], labels = 
               paste0("MDR = ", mdr), cex = 1.75)
        
        
        tmpcoor = keyALGs::plotCoor(0.5, 0.78)
        if (i == 1) Nround = 6
        empMean = keyALGs::num2str(keyALGs::Mean(X), Nround = Nround)
        text(x = tmpcoor[1], y= tmpcoor[2], labels = 
               paste0("Empirical mean = ", empMean), cex = 1.75)
        
        
        if (mdr >= mdrRangeInData[1] & mdr <= mdrRangeInData[2]) {
          lgd = c("Fitted & discretized", "Empirical") 
          fl = c("black", "black")
        }
        else  {
          lgd = c("Fitted & discretized", "Empirical, extrapolated! Fake data as surrogate!")
          fl = c("black", "red")
        }
        legend("center", legend = lgd, bty = "n", cex = 1.5, col = c("black", "gray"), 
               pch = c(NA, 15), lwd = c(1.5, NA), text.col = fl)
        
        
        Zhist = wrapAsHist(Z)
        s = (diff(range(Y$val)) / (length(Y$val) - 1)  ) / 
          (diff(range(Z$val)) / (length(Z$val) - 1) )
        l = length(Zhist$mids)
        if (Zhist$mids[l] >= 1-1e-10) 
          Zhist$density[-c(1, l)] = Zhist$density[-c(1, l)] * s
        else Zhist$density[-c(1)] = Zhist$density[-c(1)] * s
        plot(Zhist, freq = F, border = NA, col = hcol, main = "", 
             xlab = "Damage ratio", cex.lab = 2, cex.axis = 1.5, ylab = "",
             xlim = xlim, ylim = ylim, las = 1)
        lines(as.data.frame(Y), type = "h", col = "black", lwd = 1)
        
        
        zmean = keyALGs::Mean(Z)
        zmean = keyALGs::num2str(zmean, Nround = Nround)
        text(x = tmpcoor[1], y= tmpcoor[2], labels = 
               paste0("Empirical mean = ", zmean), cex = 1.5)
        legend("center", legend = "Bias corrected empirical", pch = 15, 
               col = "gray", bty = "n", cex = 1.5)
      }
      options(scipen = 0)
      dev.off()
      
    }
    
    
    # ==========================================================================
    # Make 207 plots.
    # ==========================================================================
    if (T)
    {
      tmp = unique(as.integer(round(c(
        seq(1e-5, 1e-4, by = 1e-5), seq(1e-4, 1e-2, by = 1e-4), 
        seq(1e-2, 1 - 1e-2, by = 1e-2)) * 1e5)))
      whichEmpDistrs = match(tmp, as.integer(round(targetMDRs * 1e5)))
      empDistrsToPlot = empDistrsLists$biased[whichEmpDistrs]
      empDistrsBiasCorrectedToPlot = empDistrsLists$biasCorrected[whichEmpDistrs]
      whichPMFinTable2plot = match(tmp, as.integer(round(pmftable[1,] * 1e5)))
      tablePmf2plot = pmftableToDistrList(pmftable[, whichPMFinTable2plot, drop = F])
      pdf("../figure/pmfTable-207.pdf", width = 10, height = 10 * 0.618)
      par(mar = c(4.2, 5, 0, 0), mfrow = c(2, 1), family = "serif")
      options(scipen = 999)
      for (i in 1:length(tablePmf2plot))
      {
        X = empDistrsToPlot[[i]]
        X = keyALGs::rglr(X, ngd = seq(
          X$val[1], X$val[length(X$val)], len = length(X$val)))
        Y = tablePmf2plot[[i]]
        Z = empDistrsBiasCorrectedToPlot[[i]]
        Z = keyALGs::rglr(Z, ngd = seq(
          Z$val[1], Z$val[length(Z$val)], len = length(Z$val)))
        Xhist = wrapAsHist(X)
        s = (diff(range(Y$val)) / (length(Y$val) - 1)  ) / 
          (diff(range(X$val)) / (length(X$val) - 1)  )
        l = length(Xhist$mids)
        if (Xhist$mids[l] >= 1-1e-10) 
          Xhist$density[-c(1, l)] = Xhist$density[-c(1, l)] * s
        else Xhist$density[-c(1)] = Xhist$density[-c(1)] * s
        xlim = c(0, max(c(range(Y$val), range(Xhist$breaks), range(Z$val))))
        ylim = range(c(0, max(Y$P), max(Xhist$density), max(Z$P)))
        
        
        hcol = "gray"
        plot(Xhist, freq = F, border = NA, col = hcol, main = "", 
             xlab = "", cex.lab = 2, cex.axis = 1.5, ylab = "P",
             xlim = xlim, ylim = ylim, las = 1)
        lines(as.data.frame(Y), type = "h", col = "black", lwd = 1)
        tmpcoor = keyALGs::plotCoor(0.5, 0.92)
        fullmdr = keyALGs::Mean(Y)
        Nround = 5
        
        
        mdr = keyALGs::num2str(keyALGs::Mean(Y), Nround = Nround)
        text(x = tmpcoor[1], y = tmpcoor[2], labels = 
               paste0("MDR = ", mdr), cex = 1.75)
        tmpcoor = keyALGs::plotCoor(0.5, 0.78)
        
        
        if (i == 1) Nround = 6
        empMean = keyALGs::num2str(keyALGs::Mean(X), Nround = Nround)
        text(x = tmpcoor[1], y= tmpcoor[2], labels = 
               paste0("Empirical mean = ", empMean), cex = 1.75)
        
        
        if (mdr >= mdrRangeInData[1] & mdr <= mdrRangeInData[2]) {
          lgd = c("Fitted & discretized", "Empirical") 
          fl = c("black", "black")
        }
        else  {
          lgd = c("Fitted & discretized", "Empirical, extrapolated! Fake data as surrogate!")
          fl = c("black", "red")
        }
        legend("center", legend = lgd, bty = "n", cex = 1.5, col = c("black", "gray"), 
               pch = c(NA, 15), lwd = c(1.5, NA), text.col = fl)
        
        
        Zhist = wrapAsHist(Z)
        s = (diff(range(Y$val)) / (length(Y$val) - 1)  ) / 
          (diff(range(Z$val)) / (length(Z$val) - 1) )
        l = length(Zhist$mids)
        if (Zhist$mids[l] >= 1-1e-10) 
          Zhist$density[-c(1, l)] = Zhist$density[-c(1, l)] * s
        else Zhist$density[-c(1)] = Zhist$density[-c(1)] * s
        plot(Zhist, freq = F, border = NA, col = hcol, main = "", 
             xlab = "Damage ratio", cex.lab = 2, cex.axis = 1.5, ylab = "",
             xlim = xlim, ylim = ylim, las = 1)
        lines(as.data.frame(Y), type = "h", col = "black", lwd = 1)
        
        
        zmean = keyALGs::Mean(Z)
        zmean = keyALGs::num2str(zmean, Nround = Nround)
        text(x = tmpcoor[1], y= tmpcoor[2], labels = 
               paste0("Empirical mean = ", zmean), cex = 1.5)
        legend("center", legend = "Bias corrected empirical", pch = 15, 
               col = "gray", bty = "n", cex = 1.5)
      }
      options(scipen = 0)
      dev.off()
    }
    
  }
  
  
}




# Back fit and re-discretize. Bad result, mainly due to the lowest max already
# being 1.
if (F)
{
  
  
  # ==============================================================================
  # Back fit to the pmf table.
  # ==============================================================================
  if (T)
  {
    
    
    condEmpDistrs = apply(pmftable[, -c(1, ncol(pmftable))], 2, function(x)
    {
      P = x[-(1:3)]
      list(val = seq(0, x[2], len = nrow(pmftable) - 2)[-1], P = P / sum(P))
    })
    lm1 = targetMDRs / (1 - P0)
    
    
    # startMDRtoFit = 0.05
    # startingPmf = which.min(abs(targetMDRs - startMDRtoFit))
    # cat("startingPmf =", startingPmf, "\n")
    startingPmf = 1L
    
    
    # ==========================================================================
    # Have tested sequential fitting, individual fitting. The objective function
    # values are very close.
    # ==========================================================================
    abc27 = makeInitialPoints(abcLB, abcUB)
    optRst2 = LBFGSBtrbFitList(
      # abc = abc, 
      abc = abc27, 
      lm1 = lm1,
      empDistrList = condEmpDistrs,
      abcLB = abcLB, abcUB = abcUB, 
      scaleEps = 1e-8, scaleMaxit = 100, 
      distanceFun = "likelihood",
      max_iterations = 100, maxCore = 1000, 
      RIBlib = "Numerical Recipes", 
      sequentialUpdate = startingPmf,
      # sequentialUpdate = -1, 
      hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
      epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
      max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
      ftol = 1e-4, wolfe = 0.9)
    
    
    
  }
  
  
  
  
  # ==============================================================================
  # Discretization. Round 2.
  # ==============================================================================
  if (T)
  {
    
    param = optRst$param
    d = solve_d(param, eps= 1e-8, maxit = 100)
    abcd = rbind(param[1:3, , drop = F], d)
    tailProbThreshold = 1e-3
    
    
    # ============================================================================
    # Make and lower bound the maxes.
    # ============================================================================
    maxes = pmin(1, actuar_qtrbeta(
      1 - tailProbThreshold / (1 - P0), abcd))
    maxes = pmax(0.05, maxes)
    ind = longestNonDecreasingSubseq(maxes)
    maxFun = splinefun(x = targetMDRs[ind], y = maxes[ind], method = "hyman")
    maxes = maxFun(targetMDRs)
    TrBtable = rbind(targetMDRs, maxes, P0, abcd)
    
    
    # ============================================================================
    # Discretize and tune PMF table.
    # ============================================================================
    pmftablelist = FTDFQA(
      TrBtable, supportSize = 64, 
      regridMethod = "lr",
      RIBlib = "Numerical Recipes",
      outputProbsInRows = F,
      fineDiscretizationSize = 2000,
      maxCore = 100, verbose = T, downScaleSupport = F)
    pmftable = pmftablelist$distTable
    pmftableBT = pmftablelist$distTableBeforeTune
    
    
    # ==========================================================================
    # Plot MDR vs P0, P1, maxes, cv for raw discretization.
    # ==========================================================================
    if (T)
    {
      
      # pdf("../figure/rawDiscretization-p0p1maxCv-MDR.pdf", width = 10, 
      #     height = 10 * 0.618)
      png("../figure/rawDiscretization-p0p1maxCv-MDR.png", width = 10, 
          height = 10 * 0.618, res = 120, unit = "in")
      par(mar = c(4.2, 5, 0, 0), family = "serif")
      layout(rbind(c(1, 1, 2, 2, 3, 3),
                   c(4, 5, 5, 6, 6, 7)))
      exclu = -c(1L, ncol(pmftableBT))
      # xs = 1:length(targetMDRs)
      xs = targetMDRs
      plot(x = xs, y = pmftableBT[2, exclu], type = "l", 
           col = "darkblue", bty = "L", xlim = range(xs), ylim = c(0, 1),
           ylab = "Max", cex.lab = 2, cex.axis = 1.5, las = 1, xlab = "" )
      plot(x = xs, y = pmftableBT[3, exclu], col = "darkblue", type = "l", 
           ylab = expression(P[0]), cex.lab = 2, cex.axis = 1.5, bty = "L", 
           las = 1, xlab = "" )
      plot(x = xs, y = pmftableBT[nrow(pmftableBT), exclu], 
           col = "darkblue", cex.lab = 2, cex.axis = 1.5, bty = "L", 
           las = 1, ylab = expression(P[max]), type = "l", xlab = "" )
      plot.new()
      ms = apply(pmftableBT[, exclu], 2, function(x)
      {
        sum(seq(0, x[2], len = nrow(pmftableBT) - 2) * x[-c(1, 2)])
      })
      plot(x = xs, y = ms, col = "darkblue", cex.lab = 2, cex.axis = 1.5, 
           bty = "L", las = 1, ylab = expression(mu), type = "l", xlab = "MDR" )
      sds = sqrt(apply(pmftableBT[, exclu], 2, function(x)
      {
        sum(seq(0, x[2], len = nrow(pmftableBT) - 2) ^ 2 * x[-c(1, 2)])
      }) - ms ^ 2)
      cv = sds / ms
      # cvScaled = (cv - min(cv)) / diff(range(cv))
      plot(x = xs, y = cv, col = "darkblue", cex.lab = 2, cex.axis = 1.5, 
           bty = "L", las = 1, ylab = expression(sigma/mu), type = "l", xlab = "" )
      # axis(side = 4, at = c(0, 0.5, 1), 
      #      labels = round(c(0, diff(range(cv)) / 2, max(cv)), 3), las = 1)
      dev.off()
    }
    
    
    # ==========================================================================
    # Plot MDR vs P0, P1, maxes, cv for tuned discretization.
    # ==========================================================================
    if (T)
    {
      # pdf("../figure/tunedDiscretization-p0p1maxCv-MDR.pdf", width = 10, 
      #     height = 10 * 0.618)
      png("../figure/tunedDiscretization-p0p1maxCv-MDR.png", width = 10, 
          height = 10 * 0.618, unit = "in", res = 120)
      par(mar = c(4.2, 5, 0, 0), family = "serif")
      layout(rbind(c(1, 1, 2, 2, 3, 3),
                   c(4, 5, 5, 6, 6, 7)))
      exclu = -c(1L, ncol(pmftable))
      # xs = 1:length(targetMDRs)
      xs = targetMDRs
      plot(x = xs, y = pmftable[2, exclu], type = "l", 
           col = "darkblue", bty = "L", xlim = range(xs), ylim = c(0, 1),
           ylab = "Max", cex.lab = 2, cex.axis = 1.5, las = 1, xlab = "" )
      plot(x = xs, y = pmftable[3, exclu], col = "darkblue", type = "l", 
           ylab = expression(P[0]), cex.lab = 2, cex.axis = 1.5, bty = "L", 
           las = 1, xlab = "" )
      plot(x = xs, y = pmftable[nrow(pmftable), exclu], 
           col = "darkblue", cex.lab = 2, cex.axis = 1.5, bty = "L", 
           las = 1, ylab = expression(P[max]), type = "l", xlab = "" )
      plot.new()
      ms = apply(pmftable[, exclu], 2, function(x)
      {
        sum(seq(0, x[2], len = nrow(pmftable) - 2) * x[-c(1, 2)])
      })
      plot(x = xs, y = ms, col = "darkblue", cex.lab = 2, cex.axis = 1.5, 
           bty = "L", las = 1, ylab = expression(mu), type = "l", xlab = "MDR" )
      sds = sqrt(apply(pmftable[, exclu], 2, function(x)
      {
        sum(seq(0, x[2], len = nrow(pmftable) - 2) ^ 2 * x[-c(1, 2)])
      }) - ms ^ 2)
      cv = sds / ms
      # cvScaled = (cv - min(cv)) / diff(range(cv))
      plot(x = xs, y = cv, col = "darkblue", cex.lab = 2, cex.axis = 1.5, 
           bty = "L", las = 1, ylab = expression(sigma/mu), type = "l", xlab = "" )
      # axis(side = 4, at = c(0, 0.5, 1), 
      #      labels = round(c(0, diff(range(cv)) / 2, max(cv)), 3), las = 1)
      dev.off()  
    }
    
    
    # ==========================================================================
    # Make 45 plots.
    # ==========================================================================
    if (T)
    {
      tmp = unique(as.integer(round(c(seq(1e-5, 1e-4, by = 1e-5), 
                                      seq(1e-4, 1e-3, by = 1e-4), 
                                      seq(1e-3, 1e-2, by = 1e-3), 
                                      seq(1e-2, 1e-1, by = 1e-2), 
                                      seq(1e-1, 1-1e-1, by = 1e-1)) * 1e5)))
      whichEmpDistrs = match(tmp, as.integer(round(targetMDRs * 1e5)))
      empDistrsToPlot = empDistrsLists$biased[whichEmpDistrs]
      empDistrsBiasCorrectedToPlot = empDistrs[whichEmpDistrs]
      whichPMFinTable2plot = match(tmp, as.integer(round(pmftable[1,] * 1e5)))
      tablePmf2plot = pmftableToDistrList(pmftable[, whichPMFinTable2plot, drop = F])
      pdf("../figure/pmfTable-45.pdf", width = 10, height = 10 * 0.618)
      par(mar = c(4.2, 5, 0, 0), mfrow = c(2, 1), family = "serif")
      options(scipen = 999)
      for (i in 1:length(tablePmf2plot))
      {
        X = empDistrsToPlot[[i]]
        X = keyALGs::rglr(X, ngd = seq(
          X$val[1], X$val[length(X$val)], len = length(X$val)))
        Y = tablePmf2plot[[i]]
        Z = empDistrsBiasCorrectedToPlot[[i]]
        Z = keyALGs::rglr(Z, ngd = seq(
          Z$val[1], Z$val[length(Z$val)], len = length(Z$val)))
        Xhist = wrapAsHist(X)
        s = (diff(range(Y$val)) / (length(Y$val) - 1)  ) / 
          (diff(range(X$val)) / (length(X$val) - 1)  )
        l = length(Xhist$mids)
        if (Xhist$mids[l] >= 1-1e-10) 
          Xhist$density[-c(1, l)] = Xhist$density[-c(1, l)] * s
        else Xhist$density[-c(1)] = Xhist$density[-c(1)] * s
        xlim = c(0, max(c(range(Y$val), range(Xhist$breaks), range(Z$val))))
        ylim = range(c(0, max(Y$P), max(Xhist$density), max(Z$P)))
        
        
        hcol = "gray"
        plot(Xhist, freq = F, border = NA, col = hcol, main = "", 
             xlab = "", cex.lab = 2, cex.axis = 1.5, ylab = "P",
             xlim = xlim, ylim = ylim, las = 1)
        lines(as.data.frame(Y), type = "h", col = "black", lwd = 1)
        tmpcoor = keyALGs::plotCoor(0.5, 0.92)
        fullmdr = keyALGs::Mean(Y)
        Nround = 5
        
        
        mdr = keyALGs::num2str(keyALGs::Mean(Y), Nround = Nround)
        text(x = tmpcoor[1], y = tmpcoor[2], labels = 
               paste0("MDR = ", mdr), cex = 1.75)
        
        
        tmpcoor = keyALGs::plotCoor(0.5, 0.78)
        if (i == 1) Nround = 6
        empMean = keyALGs::num2str(keyALGs::Mean(X), Nround = Nround)
        text(x = tmpcoor[1], y= tmpcoor[2], labels = 
               paste0("Empirical mean = ", empMean), cex = 1.75)
        
        
        if (mdr >= mdrRangeInData[1] & mdr <= mdrRangeInData[2])
          lgd = c("Fitted & discretized", "Empirical")
        else lgd = c("Fitted & discretized", "Empirical, extrapolated")
        legend("center", legend = lgd, bty = "n", cex = 1.5, 
               col = c("black", "gray"), pch = c(NA, 15), lwd = c(1.5, NA))
        
        
        Zhist = wrapAsHist(Z)
        s = (diff(range(Y$val)) / (length(Y$val) - 1)  ) / 
          (diff(range(Z$val)) / (length(Z$val) - 1) )
        l = length(Zhist$mids)
        if (Zhist$mids[l] >= 1-1e-10) 
          Zhist$density[-c(1, l)] = Zhist$density[-c(1, l)] * s
        else Zhist$density[-c(1)] = Zhist$density[-c(1)] * s
        plot(Zhist, freq = F, border = NA, col = hcol, main = "", 
             xlab = "Damage ratio", cex.lab = 2, cex.axis = 1.5, ylab = "",
             xlim = xlim, ylim = ylim, las = 1)
        lines(as.data.frame(Y), type = "h", col = "black", lwd = 1)
        
        
        zmean = keyALGs::Mean(Z)
        zmean = keyALGs::num2str(zmean, Nround = Nround)
        text(x = tmpcoor[1], y= tmpcoor[2], labels = 
               paste0("Empirical mean = ", zmean), cex = 1.5)
        legend("center", legend = "Bias corrected empirical", pch = 15, 
               col = "gray", bty = "n", cex = 1.5)
      }
      options(scipen = 0)
      dev.off()
      
    }
    
    
    # ==========================================================================
    # Make 207 plots.
    # ==========================================================================
    if (T)
    {
      tmp = unique(as.integer(round(c(
        seq(1e-5, 1e-4, by = 1e-5), seq(1e-4, 1e-2, by = 1e-4), 
        seq(1e-2, 1 - 1e-2, by = 1e-2)) * 1e5)))
      whichEmpDistrs = match(tmp, as.integer(round(targetMDRs * 1e5)))
      # empDistrsToPlot = empDistrs[whichEmpDistrs]
      # empDistrsBiasCorrectedToPlot = empDistrsBiasCorrected[whichEmpDistrs]
      empDistrsToPlot = empDistrsLists$biased[whichEmpDistrs]
      empDistrsBiasCorrectedToPlot = empDistrs[whichEmpDistrs]
      whichPMFinTable2plot = match(tmp, as.integer(round(pmftable[1,] * 1e5)))
      tablePmf2plot = pmftableToDistrList(pmftable[, whichPMFinTable2plot, drop = F])
      pdf("../figure/pmfTable-207.pdf", width = 10, height = 10 * 0.618)
      par(mar = c(4.2, 5, 0, 0), mfrow = c(2, 1), family = "serif")
      options(scipen = 999)
      for (i in 1:length(tablePmf2plot))
      {
        X = empDistrsToPlot[[i]]
        X = keyALGs::rglr(X, ngd = seq(
          X$val[1], X$val[length(X$val)], len = length(X$val)))
        Y = tablePmf2plot[[i]]
        Z = empDistrsBiasCorrectedToPlot[[i]]
        Z = keyALGs::rglr(Z, ngd = seq(
          Z$val[1], Z$val[length(Z$val)], len = length(Z$val)))
        Xhist = wrapAsHist(X)
        s = (diff(range(Y$val)) / (length(Y$val) - 1)  ) / 
          (diff(range(X$val)) / (length(X$val) - 1)  )
        l = length(Xhist$mids)
        if (Xhist$mids[l] >= 1-1e-10) 
          Xhist$density[-c(1, l)] = Xhist$density[-c(1, l)] * s
        else Xhist$density[-c(1)] = Xhist$density[-c(1)] * s
        xlim = c(0, max(c(range(Y$val), range(Xhist$breaks), range(Z$val))))
        ylim = range(c(0, max(Y$P), max(Xhist$density), max(Z$P)))
        
        
        hcol = "gray"
        plot(Xhist, freq = F, border = NA, col = hcol, main = "", 
             xlab = "", cex.lab = 2, cex.axis = 1.5, ylab = "P",
             xlim = xlim, ylim = ylim, las = 1)
        lines(as.data.frame(Y), type = "h", col = "black", lwd = 1)
        tmpcoor = keyALGs::plotCoor(0.5, 0.92)
        fullmdr = keyALGs::Mean(Y)
        Nround = 5
        
        
        mdr = keyALGs::num2str(keyALGs::Mean(Y), Nround = Nround)
        text(x = tmpcoor[1], y = tmpcoor[2], labels = 
               paste0("MDR = ", mdr), cex = 1.75)
        tmpcoor = keyALGs::plotCoor(0.5, 0.78)
        
        
        if (i == 1) Nround = 6
        empMean = keyALGs::num2str(keyALGs::Mean(X), Nround = Nround)
        text(x = tmpcoor[1], y= tmpcoor[2], labels = 
               paste0("Empirical mean = ", empMean), cex = 1.75)
        
        
        if (mdr >= mdrRangeInData[1] & mdr <= mdrRangeInData[2])
          lgd = c("Fitted & discretized", "Empirical")
        else lgd = c("Fitted & discretized", "Empirical, extrapolated")
        legend("center", legend = lgd, bty= "n", cex = 1.5, 
               col = c("black", "gray"), pch = c(NA, 15),
               lwd = c(1.5, NA))
        
        
        Zhist = wrapAsHist(Z)
        s = (diff(range(Y$val)) / (length(Y$val) - 1)  ) / 
          (diff(range(Z$val)) / (length(Z$val) - 1) )
        l = length(Zhist$mids)
        if (Zhist$mids[l] >= 1-1e-10) 
          Zhist$density[-c(1, l)] = Zhist$density[-c(1, l)] * s
        else Zhist$density[-c(1)] = Zhist$density[-c(1)] * s
        plot(Zhist, freq = F, border = NA, col = hcol, main = "", 
             xlab = "Damage ratio", cex.lab = 2, cex.axis = 1.5, ylab = "",
             xlim = xlim, ylim = ylim, las = 1)
        lines(as.data.frame(Y), type = "h", col = "black", lwd = 1)
        
        
        zmean = keyALGs::Mean(Z)
        zmean = keyALGs::num2str(zmean, Nround = Nround)
        text(x = tmpcoor[1], y= tmpcoor[2], labels = 
               paste0("Empirical mean = ", zmean), cex = 1.5)
        legend("center", legend = "Bias corrected empirical", pch = 15, 
               col = "gray", bty = "n", cex = 1.5)
      }
      options(scipen = 0)
      dev.off()
    }
    
    
  }
  
  
}




# ==============================================================================
# Compute ignorance score
# ==============================================================================
if (T)
{
  
  distTable = pmftablelist$distTable
  dataGrouped = aggregate(list(cdr = datResv$CDR), list(MDRind = as.integer(
    round(datResv$MDR * 1e5))), function(x) x )
  
  
  if (F) # Check P0 and P > 0
  {
    dataGrouped$cdr = lapply(dataGrouped$cdr, function(x) x[x > 1e-10])
    # dataGrouped$cdr = lapply(dataGrouped$cdr, function(x) x[x < 1e-10])
  }
  
  
  whichDistrs = match(dataGrouped$MDRind, as.integer(round(distTable[1,] * 1e5)))
  if (F) # Check if reassigning P0 will help.
  {
    distTable[3,] = c(1, P0, 0)
  }
  if (F) # Check if assigning old table's P0 and maxes will help.
  {
    distTable[2,] = oldDistTable[2, ] # Assign maxes
    distTable[3,] = oldDistTable[3, ] # Assign P0
  }
  
  
  relavantDistrs = apply(distTable[, whichDistrs, drop = F], 2, function(x)
  {
    list(val = seq(0, x[2], len = nrow(distTable) - 2L), 
         P = x[-c(1, 2)])
  })
  
  
  tscore = sum(mapply(function(x, X)
  {
    sum(igscore(X, x, redistributePwithin = T))
  }, dataGrouped$cdr, relavantDistrs, SIMPLIFY = T))
  # round(tscore, 3)
  # coverage A: -1.752 overall. Use the bias-uncorrected approach.
  # -0.915 for DR > 0
  # -0.837 for DR == 0
  # coverage C: -4.258 overall. Use the bias-uncorrected approach.
  # -0.582 for DR > 0
  # -3.676 for DR == 0
  
  
  # Check the score of the old distributions.
  if (T)
  {
    
    # oldDistTable = data.table::setDF(data.table::fread(
    #   "../tempFiles/CovC_nPts64-debora.csv"))
    # oldDistTable = t(oldDistTable)
    # dimnames(oldDistTable)= NULL
    
    
    relavantDistrsOld = apply(
      oldDistTable[, whichDistrs, drop = F], 2, function(x)
    {
      list(val = seq(0, x[2], len = nrow(oldDistTable) - 2L), 
           P = x[-c(1, 2)])
    })
    
    
    tscoreOld = sum(mapply(function(x, X)
    {
      sum(igscore(X, x, redistributePwithin = T))
    }, dataGrouped$cdr, relavantDistrsOld, SIMPLIFY = T)) #  / nrow(datResv)
    # round(tscoreOld, 3)
    # -1.824 overall.
    # -0.847 for DR > 0.
    # -0.977 for DR == 0.
    # -4.447 overall for coverage C.
    # -0.507 for DR > 0.
    # -3.94 for DR == 0.
    
  }
  
  
}
round(tscoreOld, 3); round(tscore, 3)














