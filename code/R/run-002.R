source("R/rfuns.R")
Rcpp::sourceCpp("src/amalgamation.cpp", verbose = T, cacheDir = "../tempFiles")


# ==============================================================================
# Read data in.
# ==============================================================================
load("../data/cvgAdata.Rdata")
# A claim’s damage ratio is considered to be zero if it is <= 0.002
dat = cvgData
dat$CDR[dat$CDR <= 0.002] = 0
datResv = dat[order(dat$MDR, dat$CDR), , drop = F]
rm(cvgData); gc()
# ==============================================================================




# ==============================================================================
# Read in old distribution table.
# ==============================================================================
tohoku = data.table::setDF(data.table::fread("tests/ResCovA64_Final.csv"))[-1]
tohoku = t(tohoku)
dimnames(tohoku) = NULL
tohoku = apply(tohoku, 2, function(x)
{
  list(val = seq(0, x[2], len = length(x) - 2L), P = x[-c(1, 2)])
})
tohoku = tohoku[-c(1, length(tohoku))]




# ==============================================================================
# Model P0s
# ==============================================================================
sampleSize = as.integer(round(nrow(datResv) * (2/3)))
NsampleSets = 30L
windowSizeRatio = 0.01
slidingSpeedRatioAsWindowSize = 0.01
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
  dat = dat[sample(nrow(dat), sampleSize, replace = T), , drop = F]
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




# ==============================================================================
# Bidirectional sequential fitting, no ensemble.
# ==============================================================================
if (F)
{
  
  
  source("R/rfuns.R")
  # set.seed(456)
  datNonzeroResv = datResv[datResv$CDR > 1e-10, c("MDR", "CDR")]
  datNonzeroResv = datNonzeroResv[order(
    datNonzeroResv$MDR, datNonzeroResv$CDR), , drop = F]
  # sampleSize = nrow(datNonzeroResv)
  sampleSize = as.integer(nrow(datNonzeroResv) * 0.632)
  windowSize = max(1L, as.integer(round(sampleSize * windowSizeRatio)))
  slidingSpeed = as.integer(max(1, round(
    windowSize * slidingSpeedRatioAsWindowSize)))
  windows = computeWindows(sampleSize, windowSize, start = 1L, 
                           speed = slidingSpeed)
  
  
  commonData = list(sampleSize = sampleSize, windowSize = windowSize, 
                    speed = slidingSpeed, datNozeroResv = datNozeroResv)
  
  
  f = function(seed, commonData)
  {
    set.seed(seed)
    sampleSize = commonData$sampleSize
    windowSize = commonData$windowSize
    slidingSpeed = commonData$slidingSpeed
    datNozeroResv = commonData$datNozeroResv
    
    
    dat = datNonzeroResv[sample(nrow(datNonzeroResv), sampleSize), , drop = F]
    dat = dat[order(dat$MDR), , drop = F]
    system.time({empDistrs = estimateEmpPMFs(
      mdr = dat$MDR, cdr = dat$CDR, windows, P0 = P0,
      targetMDRs = targetMDRs,
      regridMethod = "lmm",
      correctBias = F,
      nonDegenerateOldDistrs = NULL,
      sampleWeightOnOldDistrs = 0.5,
      empDistrSupportSize = 64L,
      maxCore = 1L,
      exponentialScalingForLowerMDRs = FALSE
    )})
    
    
    
  }
  
  
  
  
  set.seed(42)
  
  
  tmp = extractMain(
    empDistrs, targetMDRs, normalizeMainPart = TRUE, removeZeroPs = FALSE)
  condEmpDistrs = tmp$conditionalDistrs
  lm1 = tmp$lm1
  
  
  startMDRtoFit = mean(dat$MDR)
  startingPmf = which.min(abs(unlist(lapply(
    condEmpDistrs, function(x) sum(x[[1]] * x[[2]]))) - startMDRtoFit))
  
  
  # abc = matrix(c(6, 5, 4))
  abcLB = c(1.01, 0.1, 0.1)
  abcUB = c(30, 30, 30)
  
  
  # Manual bidirectional fitting.
  if (F)
  {
    
    rst1 = LBFGSBtrbFitList(
      abc = abc, 
      lm1 = lm1[5000:length(condEmpDistrs)],
      empDistrList = condEmpDistrs[5000:length(condEmpDistrs)],
      abcLB = abcLB, abcUB = abcUB, 
      scaleEps = 1e-8, scaleMaxit = 100, 
      distanceFun = "likelihood",
      max_iterations = 100, maxCore = 1000, 
      RIBlib = "Numerical Recipes", 
      sequentialUpdate = 1, 
      hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
      epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
      max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
      ftol = 1e-4, wolfe = 0.9)
    
    
    rst2 = LBFGSBtrbFitList(
      abc = abc, 
      lm1 = lm1[5000:1],
      empDistrList = condEmpDistrs[5000:1],
      abcLB = abcLB, abcUB = abcUB, 
      scaleEps = 1e-8, scaleMaxit = 100, 
      distanceFun = "likelihood",
      max_iterations = 100, maxCore = 1000, 
      RIBlib = "Numerical Recipes", 
      sequentialUpdate = 1, 
      hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
      epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
      max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
      ftol = 1e-4, wolfe = 0.9)
    
    
    plot(rst1$param[1,], type = "l")
    plot(rst2$param[1,], type = "l")
    
    
  }
  
  
  
  # Auto bidirectional fitting.
  if(F)
  {
    
    
    abc = matrix(c(4, 5, 6))
    # abc = makeInitialPoints(abcLB, abcUB)
    
    
    rst = LBFGSBtrbFitList(
      abc = abc, 
      lm1 = lm1,
      empDistrList = condEmpDistrs,
      abcLB = abcLB, abcUB = abcUB, 
      scaleEps = 1e-8, scaleMaxit = 100, 
      distanceFun = "likelihood",
      max_iterations = 100, maxCore = 1000, 
      RIBlib = "Numerical Recipes", 
      sequentialUpdate = 6000, 
      hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
      epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
      max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
      ftol = 1e-4, wolfe = 0.9)
    par(mfrow = c(4, 1), mar = c(0, 5, 0, 0))
    plot(rst$param[1,], type = "l", xlab = "", ylab = "", las = 1)
    plot(rst$param[2,], type = "l", xlab = "", ylab = "", las = 1)
    plot(rst$param[3,], type = "l", xlab = "", ylab = "", las = 1)
    plot(rst$fval[rst$fval < 100], type = "l", xlab = "", ylab = "", las = 1)
    
    
    par(mfrow =c(1, 1), mar = c(4.2, 5, 0, 0))
    k = 6
    tmp = as.data.frame(condEmpDistrs[[k]])
    tmpsp = seq(tmp$val[1], tmp$val[nrow(tmp)], len = 500)
    den = dtrbeta(tmpsp, rst$param[, k])
    den = den / sum(den) * (500 / nrow(tmp))
    plot(tmp, type = "h", bty = "n")
    lines(x=  tmpsp, y = den, col = "blue")
    
    
    
  }
  
  
  
  
  
  rst = LBFGSBtrbFitList(
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
  
  
  
  
}




# ==============================================================================
# Bidirectional sequential fitting, ensemble.
# ==============================================================================
if (T)
{
  
  
  source("R/rfuns.R")
  # set.seed(456)
  datNonzeroResv = datResv[datResv$CDR > 1e-10, c("MDR", "CDR")]
  datNonzeroResv = datNonzeroResv[order(
    datNonzeroResv$MDR, datNonzeroResv$CDR), , drop = F]
  windowSize = max(1L, as.integer(round(nrow(datNonzeroResv) * windowSizeRatio)))
  slidingSpeed = as.integer(max(1, round(
    windowSize * slidingSpeedRatioAsWindowSize)))
  
  
  abcLB = c(1.01, 0.1, 0.1)
  abcUB = c(30, 30, 30)
  abc = matrix(c(4, 5, 6))
  
  
  commonData = list(windowSize = windowSize, 
                    slidingSpeed = slidingSpeed, 
                    datNonzeroResv = datNonzeroResv,
                    P0 = P0, targetMDRs = targetMDRs, 
                    abcLB = abcLB, abcUB = abcUB, abc = abc)
  
  
  f = function(seed, commonData)
  {
    set.seed(seed)
    windowSize = commonData$windowSize
    slidingSpeed = commonData$slidingSpeed
    datNonzeroResv = commonData$datNonzeroResv
    targetMDRs = commonData$targetMDRs
    P0 = commonData$P0
    abcLB = commonData$abcLB
    abcUB = commonData$abcUB
    abc = commonData$abc
    
    
    dat = datNonzeroResv
    dat = dat[sample(nrow(dat), replace = F), , drop = F]
    dat = dat[order(dat$MDR), , drop = F]
    
    
    windows = computeWindows(nrow(dat), windowSize, start = 1L, 
                             speed = slidingSpeed)
    
    
    empDistrs = estimateEmpPMFs(
      mdr = dat$MDR, cdr = dat$CDR, windows, P0 = P0,
      targetMDRs = targetMDRs,
      regridMethod = "lmm",
      correctBias = T,
      nonDegenerateOldDistrs = NULL,
      sampleWeightOnOldDistrs = 0.5,
      empDistrSupportSize = 64L,
      maxCore = 1L,
      exponentialScalingForLowerMDRs = FALSE
    )
    
    
    tmp = extractMain(
      empDistrs, targetMDRs, normalizeMainPart = TRUE, removeZeroPs = FALSE)
    condEmpDistrs = tmp$conditionalDistrs
    lm1 = tmp$lm1
    
    
    startMDRtoFit = median(dat$MDR)
    startingPmf = which.min(abs(unlist(lapply(
      condEmpDistrs, function(x) sum(x[[1]] * x[[2]]))) - startMDRtoFit))
    cat("startingPmf = ", startingPmf, "\n")
    
    
    # abc = matrix(c(4, 5, 6))
    # abc = makeInitialPoints(abcLB, abcUB, 3)
    
    
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
    
    
    # ==========================================================================
    # Checking of the distance function works correctly.
    # ==========================================================================
    # ds = distances(distlist = condEmpDistrs, 
    #           param = optRst$param, scaleEps = 1e-8, 
    #           scaleMaxit = 100, 
    #           maxCore = 1000,
    #           distanceFun = "likelihood",
    #           RIBlib = "Numerical Recipes")
    # optRst$recomputedDs = ds
    
    
    optRst
  }
  
  
  source("R/CharlieMPlib.R")
  
  
  rseeds = as.list(seq(42, by = 1L, len = 1))
  ensembleRst = CharliePara(rseeds, commonData, fun = f, maxNprocess = 15)
  # ensembleRst = CharlieParaOnCluster(
    # rseeds, commonData, fun = f, maxNprocess = 288, wait = TRUE,
    # clusterHeadnodeAddress = "rscgrid163.air-worldwide.com")
  
  
  ensembleMeanEst = t(as.data.frame(lapply(1:4, function(i)
  {
    rst = apply(t(as.data.frame(lapply(
      ensembleRst, function(x) x$param[i,]))), 2, function(x)
      {
        if (i < 4) mean(x[x <= abcUB[i] & x >= abcLB[i]])
        else mean(x)
      })
    rst
  }))); dimnames(ensembleMeanEst) = NULL
  
  
  ensembleSdEst = t(as.data.frame(lapply(1:4, function(i)
  {
    rst = apply(t(as.data.frame(lapply(
      ensembleRst, function(x) x$param[i,]))), 2, function(x)
      {
        if (i < 4) sd(x[x <= abcUB[i] & x >= abcLB[i]])
        else 0
      })
    rst
  }))); dimnames(ensembleSdEst) = NULL
  
  
  # Compute the objective on the full data.
  if (F)
  {
    set.seed(42)
    
    
    dat = datNonzeroResv[sample(nrow(datNonzeroResv)), , drop = F]
    dat = dat[order(dat$MDR), , drop = F]
    
    
    windows = computeWindows(
      nrow(dat), round(nrow(dat) * windowSizeRatio), start = 1L, 
      speed = round(nrow(dat) * windowSizeRatio * slidingSpeedRatioAsWindowSize))
    
    
    empDistrs = estimateEmpPMFs(
      mdr = dat$MDR, cdr = dat$CDR, windows, P0 = P0,
      targetMDRs = targetMDRs,
      regridMethod = "lmm",
      correctBias = F,
      nonDegenerateOldDistrs = NULL,
      sampleWeightOnOldDistrs = 0.5,
      empDistrSupportSize = 64L,
      maxCore = 100L,
      exponentialScalingForLowerMDRs = FALSE
    )
    
    
    empDistrsBiasCorrected = estimateEmpPMFs(
      mdr = dat$MDR, cdr = dat$CDR, windows, P0 = P0,
      targetMDRs = targetMDRs,
      regridMethod = "lmm",
      correctBias = T,
      nonDegenerateOldDistrs = NULL,
      sampleWeightOnOldDistrs = 0.5,
      empDistrSupportSize = 64L,
      maxCore = 100L,
      exponentialScalingForLowerMDRs = FALSE
    )
    
    
    tmp = extractMain(
      empDistrs, targetMDRs, normalizeMainPart = TRUE, removeZeroPs = FALSE)
    condEmpDistrs = tmp$conditionalDistrs
    lm1 = tmp$lm1
    
    
    ds = distances(distlist = condEmpDistrs, 
                   param = ensembleMeanEst, scaleEps = 1e-8, 
                   scaleMaxit = 100, 
                   maxCore = 1000,
                   distanceFun = "likelihood",
                   RIBlib = "Numerical Recipes")
    plot(ds, type = "l")
    
    
  }
  
  
  # Make 200 plots: fitted Trb against empirical, 
  if (F)
  {
    # From 0 to 0.01, make 100 plots.
    # From 0.01 to 1, make 100 plots.
    # 1:10, seq(10L, 100L, by = 10L), seq(100L, 1000L, by = 100L)
    tmp = unique(as.integer(round(c(seq(1e-5, 1e-4, by = 1e-5), 
                              seq(1e-4, 1e-2, by = 1e-4), 
                              seq(1e-2, 1 - 1e-2, by = 1e-2)) * 1e5)))
    tmp = match(tmp, as.integer(round(targetMDRs * 1e5)))
    tmp = tmp[!is.na(tmp)]
    whichMDRtoPlot = tmp
    
    
    d = solve_d(ensembleMeanEst, eps = 1e-8, maxit = 100)
    abcd = rbind(ensembleMeanEst[1:3, ], d)
    abcdForPlot = abcd[, whichMDRtoPlot, drop = F]
    dimnames(abcdForPlot) = NULL
    
    
    empDistrsToPlot = empDistrs[whichMDRtoPlot]
    theoDistrsToPlot = mapply(function(abcd, distr)
    {
      rawDiscretization(distr[[1]][-1], abcd, lm = 1)
      # makePDFplotData(distr[[1]][2], distr[[1]][length(distr[[1]])], 
      #                 abcd, lm = 1, npoint = 500, normalize = TRUE)
    }, as.data.frame(abcdForPlot), empDistrsToPlot, SIMPLIFY = F)
    
    
    # range(unlist(lapply(theoDistrsToPlot, function(x) x[[1]][1])) - 
    #         unlist(lapply(empDistrsToPlot, function(x) x[[1]][2])))
    
    
    pdf("../figure/ensemble-emp-theo-rawDiscretization.pdf", 
        width = 10, height = 10 * (9/16))
    par(mar = c(4.2, 5, 0, 0), mfrow = c(1, 1), family = "serif")
    for (i in 1:length(empDistrsToPlot))
    {
      X = empDistrsToPlot[[i]]
      Y = theoDistrsToPlot[[i]]
      Y$P = Y$P * (1 - X$P[1])
      ylim = c(0, max(max(X$P), max(Y$P)))
      Xhist = wrapAsHist(X)
      Yscaled = as.data.frame(Y)
      Yscaled$P = Yscaled$P * (length(Y$val) / length(X[[1]]))
      # ylim = c(0, 1)
      plot(Xhist, ylim = ylim, 
           border = NA,
           col = "gray", freq = F, main = "",
           #, xlim = c(0, 0.1)
           cex.lab = 2, cex.axis = 1.5
           , xlab = "Damage ratio"
           )
      lines(Yscaled, col = "blue", lwd = 1)
    }
    dev.off()
    
    
  }
  
  
  # ============================================================================
  # Make PMF table.
  # ============================================================================
  if (T)
  {
    
    
    d = solve_d(ensembleMeanEst, eps= 1e-8, maxit = 100)
    abcd = rbind(ensembleMeanEst[1:3, , drop = F], d)
    tailProbThreshold = 1e-5
    maxes = pmin(1, actuar_qtrbeta(
      1 - tailProbThreshold / (1 - P0), abcd))
    
    
    # ==========================================================================
    # Make and lower bound the maxes.
    # ==========================================================================
    maxes = pmax(0.05, maxes)
    ind = longestNonDecreasingSubseq(maxes)
    maxFun = splinefun(x = targetMDRs[ind], y = maxes[ind], method = "hyman")
    maxes = maxFun(targetMDRs)
    TrBtable = rbind(targetMDRs, maxes, P0, abcd)
    
    
    # ==========================================================================
    # Discretize and tune PMF table.
    # ==========================================================================
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
      
      
      pdf("../figure/rawDiscretization-p0p1maxCv-MDR.pdf", width = 10, 
          height = 10 * 0.618)
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
      pdf("../figure/tunedDiscretization-p0p1maxCv-MDR.pdf", width = 10, 
          height = 10 * 0.618)
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
      empDistrsToPlot = empDistrs[whichEmpDistrs]
      empDistrsBiasCorrectedToPlot = empDistrsBiasCorrected[whichEmpDistrs]
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
        
        
        legend("center", legend = c("Fitted & discretized", "Empirical"), 
               bty= "n", cex = 1.5, col = c("black", "gray"), pch = c(NA, 15),
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
    
    
    # ==========================================================================
    # Make 207 plots.
    # ==========================================================================
    if (T)
    {
      tmp = unique(as.integer(round(c(
        seq(1e-5, 1e-4, by = 1e-5), seq(1e-4, 1e-2, by = 1e-4), 
        seq(1e-2, 1 - 1e-2, by = 1e-2)) * 1e5)))
      whichEmpDistrs = match(tmp, as.integer(round(targetMDRs * 1e5)))
      empDistrsToPlot = empDistrs[whichEmpDistrs]
      empDistrsBiasCorrectedToPlot = empDistrsBiasCorrected[whichEmpDistrs]
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
        
        
        legend("center", legend = c("Fitted & discretized", "Empirical"), 
               bty= "n", cex = 1.5, col = c("black", "gray"), pch = c(NA, 15),
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
  
  
  # ============================================================================
  # Back fitting to the PMF table.
  # ============================================================================
  if (T)
  {
    condEmpDistrs = apply(pmftable[, -c(1, ncol(pmftable))], 2, function(x)
    {
      P = x[-(1:3)]
      list(val = seq(0, x[2], len = nrow(pmftable) - 2)[-1], P = P / sum(P))
    })
    lm1 = targetMDRs / (1 - P0)
    
    
    startMDRtoFit = 0.05
    startingPmf = which.min(abs(targetMDRs - startMDRtoFit))
    cat("startingPmf = ", startingPmf, "\n")
    
    
    # ==========================================================================
    # Have tested sequential fitting, individual fitting. The objective function
    # values are very close.
    # ==========================================================================
    # abc27 = makeInitialPoints(abcLB, abcUB)
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
      # sequentialUpdate = -1, 
      hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
      epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
      max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
      ftol = 1e-4, wolfe = 0.9)
    
    
    
    
  }
  
  
  
  
}

