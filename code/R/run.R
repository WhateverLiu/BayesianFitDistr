source("R/rfuns.R")
Rcpp::sourceCpp("src/amalgamation.cpp", verbose = T)


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


# modelP0 = function(mdr, cdr, windows, targetMDRs, 
#                    eps = 1e-10, P0modeltype = "hyman")
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
    eps = 1e-10, P0modeltype = "hyman")
}
P0 = rowMeans(as.data.frame(p0models))
plot(x = targetMDRs, y = P0, type = "l")




source("R/rfuns.R")
set.seed(456)
datNonzeroResv = datResv[datResv$CDR > 1e-10, c("MDR", "CDR")]
datNonzeroResv = datNonzeroResv[order(
  datNonzeroResv$MDR, datNonzeroResv$CDR), , drop = F]
sampleSize = as.integer(round(nrow(datNonzeroResv) * (2/3)))
windowSize = max(1L, as.integer(round(sampleSize * windowSizeRatio)))
slidingSpeed = as.integer(max(1, round(
  windowSize * slidingSpeedRatioAsWindowSize)))
windows = computeWindows(sampleSize, windowSize, start = 1L, speed = slidingSpeed)


sequentialUpdate = T
startMDRtoFit = mean(datResv$MDR)
distanceFun = "likelihood"
NsampleSets = 30L
ensembleRst = list()
multipleInits = T
for (i in 1:NsampleSets)
{
  cat("Sample set index =", i, "\n")
  
  
  dat = datNonzeroResv
  ind = sample(nrow(dat), sampleSize, replace = F)
  
  
  dat = dat[ind, , drop = F]
  dat = dat[order(dat$MDR), , drop = F]
  
  
  # system.time({p0andEmpDistrs = modelP0andMakeEmpiricalDistributions(
  system.time({p0andEmpDistrs = estimateEmpPMFs(
    mdr = dat$MDR, cdr = dat$CDR, windows, P0 = P0,
    targetMDRs = targetMDRs,
    regridMethod = "lmm",
    correctBias = FALSE,
    # correctBias = TRUE,
    nonDegenerateOldDistrs = NULL,
    # nonDegenerateOldDistrs = tohoku,
    sampleWeightOnOldDistrs = 0.99,
    empDistrSupportSize = 64L,
    maxCore = 1000L,
    exponentialScalingForLowerMDRs = FALSE
  )})
  
  
  # Plot all empirical PMFs.
  if (F)
  {
    
    mdrInterest = unique(
      c(1:10, seq(10L, 100L, by = 10), seq(100L, 1000L, by = 100L),
        seq(1000L, 10000L, by = 1000L), seq(10000L, 100000L, by = 10000L))) / 1e5
    mdrInterest = mdrInterest[mdrInterest < 1]
    ind = match(as.integer(round(mdrInterest * 1e6)), 
                as.integer(round(p0andEmpDistrs$targetMDRs * 1e6)))
    tmp = p0andEmpDistrs$empDistrs[ind]
    ylim = c(0, max(unlist(lapply(tmp, function(x) 
      max(x[[2]][-c(1,length(x[[2]]))])))))
    
    
    pdf("../figure/tmpempDistrs.pdf", width = 13, height = 13 * (9/16))
    options(scipen = 999)
    par(mar = c(2, 3.3, 0, 0), mfrow = c(5, 9), family = "serif")
    for (i in 1:length(tmp))
    {
      if (i == 1) yaxt = 's' else yaxt = 'n'
      if (i == 37) xaxt = 's' else xaxt = 'n'
      h = wrapAsHist(tmp[[i]])
      plot(h, main = "",
           # type = "h", bty = "L", 
           xlim = c(0, 1), 
           ylim = ylim, cex.lab = 1.5, cex.axis = 1.5, xlab = "", col = "gray",
           ylab = "", las = 1, xaxt = xaxt, yaxt = yaxt, freq = F, border = NA)
      if (i == 1) title(xlab = "DR", line = 0.75, cex.lab = 1.5)
      if (i == 10) title(ylab = "P", line = 1, cex.lab = 1.5)
      if (i == 1) lgd = paste0("MDR = \n", mdrInterest[i])
      else lgd = paste0(mdrInterest[i])
      legend("top", legend = lgd, bty= "n", cex = 1.5)
    }
    options(scipen = 0)
    dev.off()
    
    
    pdf("../figure/tmpempDistrsMainIndividualScope.pdf", width = 13, height = 13 * (9/16))
    options(scipen = 999)
    par(mar = c(2, 3.3, 0, 0), mfrow = c(5, 9), family = "serif")
    for (i in 1:length(tmp))
    {
      if (i == 1) yaxt = 's' else yaxt = 'n'
      if (i == 37) xaxt = 's' else xaxt = 'n'
      h = wrapAsHist(tmp[[i]])
      ylim = c(0, max(h$density[-1]))
      plot(h, main = "",
           # type = "h", bty = "L", 
           # xlim = c(0, 1), 
           ylim = ylim,
           cex.lab = 1.5, cex.axis = 1.5, xlab = "", col = "gray",
           ylab = "", las = 1, xaxt = xaxt, yaxt = yaxt, freq = F, border = NA)
      if (i == 1) title(xlab = "DR", line = 0.75, cex.lab = 1.5)
      if (i == 10) title(ylab = "P", line = 1, cex.lab = 1.5)
      if (i == 1) lgd = paste0("MDR = \n", mdrInterest[i])
      else lgd = paste0(mdrInterest[i])
      legend("top", legend = lgd, bty= "n", cex = 1.5)
    }
    options(scipen = 0)
    dev.off()
    
    
    
  }
  
  
  # ==============================================================================
  # Create initial parameters.
  # ==============================================================================
  ndiv = 3L
  abcLB = c(1.01, 0.1, 0.1)
  abcUB = c(30, 30, 30)
  
  
  if (multipleInits)
  {
    abc = apply(rbind(abcLB, abcUB), 2, function(x)
    {
      exp(seq(log(x[1]), log(x[2]), len = ndiv + 2L)[2L:(2L + ndiv - 1L)])
    })
    abc = t(expand.grid(a = abc[,1], b = abc[,2], c = abc[,3]))
  }
  else 
  {
    abc = as.matrix(c(4, 5, 6))
  }
  
  
  tmp = extractMain(p0andEmpDistrs, targetMDRs, 
                    normalizeMainPart = TRUE, removeZeroPs = FALSE)
  condEmpDistrs = tmp$conditionalDistrs
  lm1 = tmp$lm1
  
  
  startingPmf = which.min(abs(unlist(lapply(
    p0andEmpDistrs, function(x) sum(x[[1]] * x[[2]]))) - startMDRtoFit))
  
  
  # Optimization.
  if (T)
  {
    if (!sequentialUpdate) startingPmf = 1L
    
    
    secondHalf = startingPmf:length(condEmpDistrs)
    system.time({optRstSecondHalf = LBFGSBtrbFitList(
      abc = abc, 
      lm1 = lm1[secondHalf],
      empDistrList = condEmpDistrs[secondHalf], 
      abcLB = abcLB, abcUB = abcUB, 
      scaleEps = 1e-8, scaleMaxit = 100, 
      distanceFun = distanceFun,
      max_iterations = 100, maxCore = 1000, 
      RIBlib = "Numerical Recipes", 
      sequentialUpdate = sequentialUpdate, 
      hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
      epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
      max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
      ftol = 1e-4, wolfe = 0.9)
    })
    
    
    if (startingPmf != 1L)
    {
      firstHalf = startingPmf:1  
      system.time({optRstFirstHalf = LBFGSBtrbFitList(
        abc = abc, 
        lm1 = lm1[firstHalf],
        empDistrList = condEmpDistrs[firstHalf], 
        abcLB = abcLB, abcUB = abcUB, 
        scaleEps = 1e-8, scaleMaxit = 100, 
        distanceFun = distanceFun,
        max_iterations = 100, maxCore = 1000, 
        RIBlib = "Numerical Recipes", 
        sequentialUpdate = sequentialUpdate, 
        hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
        epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
        max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
        ftol = 1e-4, wolfe = 0.9)
      })
      optRst = list()
      optRst$param = cbind(optRstFirstHalf$param[, firstHalf], 
                           optRstSecondHalf$param[, -1])
      optRst$fval = c(optRstFirstHalf$fval[firstHalf], 
                      optRstSecondHalf$fval[-1])
      optRst$niter = c(optRstFirstHalf$niter[firstHalf], 
                       optRstSecondHalf$niter[-1])
    }
    else optRst = optRstSecondHalf

    
  }
  
  
  ensembleRst[[length(ensembleRst) + 1]] = list(
    optRst = optRst, empMeans = unlist(
      lapply(p0andEmpDistrs, function(x) keyALGs::Mean(x))))
}


# save(ensembleRst, file = "tmp/ensembleRstSeq.Rdata")
save(ensembleRst, file = "tmp/ensembleRstNonSeqMutiInit.Rdata")


load("tmp/ensembleRstSeq.Rdata")
optRstNonSeq = new.env()
load("tmp/ensembleRstNonSeqMutiInit.Rdata", envir = optRstNonSeq)


# Plot MDR vs mean and conditional mean.
if (T)
{
  
  
  pdf("../figure/P0-empirical-mean-MDR.pdf", width = 10, height = 10 * 0.5)
  # png("../figure/P0-empirical-mean-MDR.png", width = 10, height = 10 * 0.5,
      # res = 120, units = "in")
  par(mfrow = c(1,2), family = "serif", mar = c(4.2, 5, 0, 0))
  plot(x = targetMDRs, y = P0, type ="l",
    bty = "L", xlab = "MDR", ylab = expression(P[0]), 
    las = 1, cex.lab = 2, cex.axis =1.5, lwd = 1)
  lines(0, col = scales::alpha("white", 0))
  rect(xleft = mdrRangeInData[1], ybottom = -1, xright = mdrRangeInData[2], 
       ytop = 10, col = scales::alpha("gray", 0.3), border = NA)
  legend("topright", legend = paste0(
    "MDR range\nin data =\n", 
    paste0('[', paste0(round(mdrRangeInData, 4), collapse = ", "), ']')), 
    pch = 15, col = "gray", cex = 1.5, bty = "n")
  m = rowMeans(as.data.frame(lapply(ensembleRst, function(x) x$empMeans)))
  plot(x = targetMDRs, y = m, 
    type = "l", col = "blue", xlab = "", ylab = "Empirical Mean", las = 1,
    cex.lab = 2, cex.axis =1.5, bty = "L", lwd = 1, xlim = c(0, 1), ylim = c(0, 1))
  lines(x = c(0, 1), y = c(0, 1), col = "red", lty = 2)
  rect(xleft = mdrRangeInData[1], ybottom = -1, xright = mdrRangeInData[2], 
       ytop = 10, col = scales::alpha("gray", 0.3), border = NA )
  dev.off()
  
  
}


# Compare non-sequential fitting vs sequential fitting
if (T)
{
  
  # Metropolis Hasting (Simulated Annealing)
  # Non sequential Quasi-Newton without ensemble estimation, 
  #   27 initializations, select the best.
  # Non sequential Quasi-Newton with bagging, take the average. 
  #   27 initializations, select the best outcome.
  # Sequential Quasi-Newton fitting without bagging.
  # Sequential fitting with bagging.
  png("../figure/abcFromDifferenceOpt.png", width = 16, height = 16 * 0.5,
      res = 360, units = "in")
  # par(mar = c(4.2, 5, 0, 0), mfrow = c(6, 4))
  # par(mar = c(4.2, 5, 0, 0), family = "serif")
  par(family = "serif")
  tmp = t(matrix(c(rep(1, 4), 2:5, 6:9, 10:13, rep(14, 4),
    15:18, 19:22, 23:26), nrow = 4))
  layout(mat = tmp, heights = c(1/2, 1, 1, 1, 1/2, 1, 1, 1), widths = c(1, 1, 1, 1))
  # plot.new(); plot.new(); plot.new(); plot.new()
  par(mar = c(0, 0, 0, 0))
  plot(0, bty = "n", xaxt = "n", yaxt = "n", 
       xlab = "", ylab = "", col = "white")
  legend("center", legend = "Fitting individually", cex = 2, bty = "n")
  
  
  par(mar = c(4.2, 5, 0, 0))
  tmpf = function(X, i)
  {
    X = t(as.data.frame(X))
    dimnames(X) = NULL
    apply(X, 2, function(x) mean(x[x >= abcLB[i] & x <= abcUB[i]]))
  }
  
  
  tmpf2 = function(X)
  {
    X = t(as.data.frame(X))
    dimnames(X) = NULL
    apply(X, 2, function(x) mean( x[x > -100 & x < 100] ))
  }
  
  
  abc = optRstNonSeq$ensembleRst[[1]]$optRst$param
  u = abc[1,]; u = u[u >= abcLB[1] & u <= abcUB[1]]
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
       xlab = "", ylab = "Full data", las = 1, cex.lab = 1.5
       , ylim = c(abcLB[1], abcUB[1]) 
       ,col = "gray"
       )
  legend("top", legend = "a", bty = "n", cex = 2)
  u = abc[2,]; u = u[u >= abcLB[2] & u <= abcUB[2]]
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
       xlab = "", ylab = "", las = 1
       , ylim = c(abcLB[2], abcUB[2])
       ,col = "gray")
  legend("top", legend = "b", bty = "n", cex = 2)
  u = abc[3,]; u = u[u >= abcLB[3] & u <= abcUB[3]]
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n", 
       xlab = "", ylab = "", las = 1
       , ylim = c(abcLB[3], abcUB[3])
       ,col = "gray")
  legend("top", legend = "c", bty = "n", cex = 2)
  u = optRstNonSeq$ensembleRst[[1]]$optRst$fval
  u = u[u > -100 & u < 100]
  plot(optRstNonSeq$ensembleRst[[1]]$optRst$fval, 
       type = "l", bty = "L", cex.axis = 1.5, xaxt = "n", 
       xlab = "", ylab = "", las = 1, ylim = c(-0.12, 0.12))
  tmp = keyALGs::plotCoor(0.5, 0.9)
  text(x = tmp[1], y= tmp[2], 
       labels = "Negative log-likelihood", cex = 1.5)
  
  
  set.seed(1234)
  Ndist = ncol(optRstNonSeq$ensembleRst[[1]]$optRst)
  cols = scales::alpha(randomcoloR::distinctColorPalette(
    length(optRstNonSeq$ensembleRst)), 0.35)
  
  
  getParam = function(rst, whichParm = 1) # 4 gives fval.
  {
    k = whichParm
    if (k < 4)
    {
      lapply(rst, function(x) 
      {
        u = x$optRst$param[k, ]
        ind = which(u >= abcLB[k] & u <= abcUB[k])
        data.frame(ind, u[ind])
      })  
    }
    else 
    {
      lapply(rst, function(x) 
      {
        u = x$optRst$fval
        ind = which(u > -100 & u < 100)
        data.frame(ind, u[ind])
      })
    }
    
  }
  
  
  # Plot the 30 lines for a, b, c in non-seq opt result.
  if (T)
  {
    for (k in 1:4)
    {
      tmp = getParam(optRstNonSeq$ensembleRst, k)
      for (j in 1:length(tmp))
      {
        u = tmp[[j]]
        if (j == 1) plt = plot else plt = lines
        if (j == 1 & k == 1) ylab = "Ensemble" else ylab = ""
        if (k != 4) ylim = c(0, 30) else ylim = c(-0.12, 0.12)
        plt(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
            xlab = "", ylab = ylab, las = 1, col = cols[j], 
            cex.lab = 1.5, ylim = ylim)
      }  
    }
  }
  
  
  tmp = optRstNonSeq$ensembleRst
  u = tmpf(lapply(tmp, function(x) x$optRst$param[1,]), 1)
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
       xlab = "", ylab = "Average", las = 1, cex.lab = 1.5
       , ylim = c(abcLB[1], abcUB[1])
       # ,col = "gray"
       )
  u = tmpf(lapply(tmp, function(x) x$optRst$param[2,]), 2)
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
       xlab = "", ylab = "", las = 1
       , ylim = c(abcLB[2], abcUB[2])
       # ,col = "gray"
       )
  u = tmpf(lapply(tmp, function(x) x$optRst$param[3,]), 3)
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
       xlab = "", ylab = "", las = 1
       , ylim = c(abcLB[3], abcUB[3])
       # ,col = "gray"
       )
  u = tmpf2( lapply(tmp, function(x) x$optRst$fval)  )
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
       xlab = "", ylab = "", las = 1, ylim = c(-0.12, 0.12))
  
  
  par(mar = c(0, 0, 0, 0))
  plot(0, bty = "n", xaxt = "n", yaxt = "n", 
       xlab = "", ylab = "", col = "white")
  legend("top", legend = "Bidirectional sequential fitting", cex = 2, bty = "n")
  par(mar = c(4.2, 5, 0, 0))
  
  
  abc = ensembleRst[[1]]$optRst$param
  u = abc[1,]; u = u[u >= abcLB[1] & u <= abcUB[1]]
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
       xlab = "", las = 1, ylab = "Full data", cex.lab = 1.5
       , ylim = c(abcLB[1], abcUB[1])
       )
  u = abc[2,]; u = u[u >= abcLB[2] & u <= abcUB[2]]
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
       xlab = "", ylab = "", las = 1
       , ylim = c(abcLB[2], abcUB[2])
       )
  u = abc[3,]; u = u[u >= abcLB[3] & u <= abcUB[3]]
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n", 
       xlab = "", ylab = "", las = 1
       , ylim = c(abcLB[3], abcUB[3])
       )
  u = ensembleRst[[1]]$optRst$fval
  u = u[u > -100 & u < 100]
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n", 
       xlab = "", ylab = "", las = 1, ylim = c(-0.12, 0.12))
  
  
  # Plot the 30 lines for a, b, c in non-seq opt result.
  if (T)
  {
    for (k in 1:4)
    {
      tmp = getParam(ensembleRst, k)
      for (j in 1:length(tmp))
      {
        if (k == 1 & j == 1) ylab = "Ensemble" else ylab = ""
        u = tmp[[j]]
        if (j == 1) plt = plot else plt = lines
        if (k != 4) ylim = c(0, 30)
        else ylim = c(-0.12, 0.12)
        plt(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
            xlab = "", ylab = ylab, cex.lab = 1.5, 
            las = 1, col = cols[j], ylim = ylim)
      }  
    }
  }
  
  
  tmp = ensembleRst
  u = tmpf(lapply(tmp, function(x) x$optRst$param[1,]), 1)
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
       xlab = "", las = 1
       , ylim = c(abcLB[1], abcUB[1]), ylab = "Average", cex.lab = 1.5
       )
  axis(side = 1, at = c(1, 10000, 18999), labels = c(1, 10000, 18999), 
       cex.axis = 1.5)
  title(xlab = "MDR index", cex.lab = 2)
  u = tmpf(lapply(tmp, function(x) x$optRst$param[2,]), 2)
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
       xlab = "", ylab = "", las = 1
       , ylim = c(abcLB[2], abcUB[2])
       )
  u = tmpf(lapply(tmp, function(x) x$optRst$param[3,]), 3)
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
       xlab = "", ylab = "", las = 1
       , ylim = c(abcLB[3], abcUB[3])
       )
  u = tmpf2( lapply(tmp, function(x) x$optRst$fval)  )
  plot(u, type = "l", bty = "L", cex.axis = 1.5, xaxt = "n",
       xlab = "", ylab = "", las = 1, ylim = c(-0.12, 0.12))
  
  
  
  dev.off()
  
  
}




# ==============================================================================
# Test mixing with old distributions.
# ==============================================================================
tohoku = data.table::setDF(data.table::fread("tests/ResCovA64_Final.csv"))[-1]
tohoku = t(tohoku)
dimnames(tohoku) = NULL
tohoku = apply(tohoku, 2, function(x)
{
  list(val = seq(0, x[2], len = length(x) - 2L), P = x[-c(1, 2)])
})
tohoku = tohoku[-c(1, length(tohoku))]
tohokuP0 = unlist(lapply(tohoku, function(x) x[[2]][1]))
sampleWeights = c(0.99, 0.9, 0.8, 0.5, 0.2, 0.1, 0.01)
newP0list = lapply(sampleWeights, function(w) (1 - w) * P0 + tohokuP0 * w)


dat = datResv[c("MDR", "CDR")]
dat = dat[dat$CDR > 0, , drop = F]
set.seed(789)
dat = dat[sample(nrow(dat), nrow(dat)), , drop = F]
dat = dat[order(dat$MDR), , drop = F]


sampleSize = nrow(dat)
windowSizeRatio = 0.01
slidingSpeedRatioAsWindowSize = 0.01
windowSize = max(1L, as.integer(round(sampleSize * windowSizeRatio)))
slidingSpeed = as.integer(max(1, round(
  windowSize * slidingSpeedRatioAsWindowSize)))
windows = computeWindows(sampleSize, windowSize, 
                         start = 1L, speed = slidingSpeed)

sequentialUpdate = T
startMDRtoFit = mean(datResv$MDR)
distanceFun = "likelihood"
# NsampleSets = 30L
multipleInits = F


ensembleRst = list()
for (i in 1:length(newP0list))
{
  cat(i, "")
  
  
  newP0 = newP0list[[i]]
  system.time({p0andEmpDistrs = estimateEmpPMFs(
    mdr = dat$MDR, cdr = dat$CDR, windows, P0 = newP0,
    targetMDRs = targetMDRs,
    regridMethod = "lmm",
    correctBias = FALSE,
    # correctBias = TRUE,
    # nonDegenerateOldDistrs = NULL,
    nonDegenerateOldDistrs = tohoku,
    # sampleWeightOnOldDistrs = ,
    sampleWeightOnOldDistrs = sampleWeights[i],
    empDistrSupportSize = 64L,
    maxCore = 1000L,
    exponentialScalingForLowerMDRs = FALSE
  )})
  
  
  # ==============================================================================
  # Create initial parameters.
  # ==============================================================================
  ndiv = 3L
  abcLB = c(1.01, 0.1, 0.1)
  abcUB = c(30, 30, 30)
  
  
  if (multipleInits)
  {
    abc = apply(rbind(abcLB, abcUB), 2, function(x)
    {
      exp(seq(log(x[1]), log(x[2]), len = ndiv + 2L)[2L:(2L + ndiv - 1L)])
    })
    abc = t(expand.grid(a = abc[,1], b = abc[,2], c = abc[,3]))
  } else 
  {
    abc = as.matrix(c(4, 5, 6))
  }
  
  
  tmp = extractMain(p0andEmpDistrs, targetMDRs, 
                    normalizeMainPart = TRUE, removeZeroPs = FALSE)
  condEmpDistrs = tmp$conditionalDistrs
  lm1 = tmp$lm1
  
  
  startingPmf = which.min(abs(unlist(lapply(
    p0andEmpDistrs, function(x) sum(x[[1]] * x[[2]]))) - startMDRtoFit))
  
  
  # Optimization.
  if (T)
  {
    if (!sequentialUpdate) startingPmf = 1L
    
    
    secondHalf = startingPmf:length(condEmpDistrs)
    system.time({optRstSecondHalf = LBFGSBtrbFitList(
      abc = abc, 
      lm1 = lm1[secondHalf],
      empDistrList = condEmpDistrs[secondHalf], 
      abcLB = abcLB, abcUB = abcUB, 
      scaleEps = 1e-8, scaleMaxit = 100, 
      distanceFun = distanceFun,
      max_iterations = 100, maxCore = 1000, 
      RIBlib = "Numerical Recipes", 
      sequentialUpdate = sequentialUpdate, 
      hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
      epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
      max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
      ftol = 1e-4, wolfe = 0.9)
    })
    
    
    if (startingPmf != 1L)
    {
      firstHalf = startingPmf:1  
      system.time({optRstFirstHalf = LBFGSBtrbFitList(
        abc = abc, 
        lm1 = lm1[firstHalf],
        empDistrList = condEmpDistrs[firstHalf], 
        abcLB = abcLB, abcUB = abcUB, 
        scaleEps = 1e-8, scaleMaxit = 100, 
        distanceFun = distanceFun,
        max_iterations = 100, maxCore = 1000, 
        RIBlib = "Numerical Recipes", 
        sequentialUpdate = sequentialUpdate, 
        hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
        epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
        max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
        ftol = 1e-4, wolfe = 0.9)
      })
      optRst = list()
      optRst$param = cbind(optRstFirstHalf$param[, firstHalf], 
                           optRstSecondHalf$param[, -1])
      optRst$fval = c(optRstFirstHalf$fval[firstHalf], 
                      optRstSecondHalf$fval[-1])
      optRst$niter = c(optRstFirstHalf$niter[firstHalf], 
                       optRstSecondHalf$niter[-1])
    } else
    {
      optRst = optRstSecondHalf
    }
      
    
    
  }
  
  
  ensembleRst[[length(ensembleRst) + 1]] = list(
    optRst = optRst, empMeans = unlist(
      lapply(p0andEmpDistrs, function(x) keyALGs::Mean(x))))
  
  
  
}


# Compare PMFs from Tohoku and fitted.
if (T)
{
  
  mdrInterested = unique(as.integer(round(
    c(1:10, seq(10, 100, by = 10), seq(1e2, 1e3, by = 1e2),
    seq(1e3, 1e4, by = 1e3), seq(1e4, 1e5, by =1e4)))))
  mdrInterested = mdrInterested[-length(mdrInterested)] / 1e5
  
  
  mdrInterestedInd = match(as.integer(round(mdrInterested * 1e6)),
                       as.integer(round(targetMDRs * 1e6)))
  
  
  tohokuInterested = tohoku[mdrInterestedInd]
  tohokuMaxes = unlist(lapply(tohokuInterested, function(x) 
    x[[1]][length(x[[1]])]))
  
  
  distlist = mapply(function(x, y)
  {
    list(mdr = mdrInterested,
         maxes = tohokuMaxes,
         p0 = y[mdrInterestedInd],
         param = x$optRst$param[, mdrInterestedInd, drop = F])
  }, ensembleRst, newP0list, SIMPLIFY = F)
  
  
  tmpf = function(X)
  {
    ds = solve_d(X$param, eps = 1e-8, maxit = 100)
    pm = rbind(X$param[-4, ], ds)
    mx = X$maxes
    mdr = X$mdr
    rst = list()
    for (i in 1:length(mx))
    {
      sp = seq(0, mx[i], len = 64)
      delta = sp[2] - sp[1]
      lb = sp[-1] - delta / 2
      ub = lb + delta
      ub = ub[-length(ub)]
      p2 = c(actuar_ptrbeta(ub, pm[, i, drop = F]), 1)
      p1 = actuar_ptrbeta(lb, pm[, i, drop = F])
      pmf = p2 - p1
      pmf = pmf / sum(pmf) * (1 - X$p0[i])
      tmpP = c(X$p0[i], pmf)
      rst[[i]] = list(val = sp, P = tmpP )
      
      
      # pmf = list(val = sp, P = tmpP )
      # m = sum(pmf[[1]] * pmf[[2]])
      # pmf = upScaleAndBoundPMF(pmf, 1, mdr[i] / m, regridMethod = "lr")
      # rst[[i]] = pmf[[1]]
    }; rst
  }
  
  
  pmflist = lapply(distlist, function(X) tmpf(X))
  
  
  pdf("../figure/tohokuVSnz.pdf", width = 10, height = 10 * 0.45)
  options(scipen = 999)
  par(mar = c(0.3, 0.3, 0.3, 0.3), family = "serif", mfrow = c(5, 9))
  layout(t(matrix(c(rep(1, 9), 2:46), nrow = 9)), 
         heights = c(1, rep(2, 5)))
  for (i in length(pmflist):1)
  {
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = paste0("Tohoku weight = ", sampleWeights[i]),
           bty = "n", cex = 2)
    for (j in 1:length(pmflist[[i]]))
    {
      tok = tohokuInterested[[j]]
      h = wrapAsHist(tok)
      s = wrapAsHist(pmflist[[i]][[j]])
      ylim = range(c(0, max(s$density), max(h$density)))
      plot(h, freq = F, border = NA, col = scales::alpha("red", 0.5), 
           main = "", xaxt = "n", yaxt = "n", ylim = ylim)
      lines(s, type = "h", border = NA, col = scales::alpha("blue", 0.5), 
            freq = F)
      if (j == 1) lgd = paste0("MDR =\n", mdrInterested[j])
      else lgd = mdrInterested[j]
      legend("top", legend = lgd, cex = 1.5, bty = "n")
    }
  }
  options(scipen = 0)
  dev.off()
  
  
  pdf("../figure/tohokuVSnz-mainPart.pdf", width = 10, height = 10 * 0.45)
  options(scipen = 999)
  par(mar = c(0.3, 0.3, 0.3, 0.3), family = "serif", mfrow = c(5, 9))
  layout(t(matrix(c(rep(1, 9), 2:46), nrow = 9)), 
         heights = c(1, rep(2, 5)))
  for (i in length(pmflist):1)
  {
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = paste0("Tohoku weight = ", sampleWeights[i]),
           bty = "n", cex = 2)
    for (j in 1:length(pmflist[[i]]))
    {
      tok = tohokuInterested[[j]]
      h = wrapAsHist(tok)
      s = wrapAsHist(pmflist[[i]][[j]])
      l = -c(1, length(s$density))
      ylim = range(c(0, max(s$density[l]), max(h$density[l])))
      plot(h, freq = F, border = NA, col = scales::alpha("red", 0.5), 
           main = "", xaxt = "n", yaxt = "n", ylim = ylim)
      lines(s, type = "h", border = NA, col = scales::alpha("blue", 0.5), 
            freq = F)
      if (j == 1) lgd = paste0("MDR =\n", mdrInterested[j])
      else lgd = mdrInterested[j]
      legend("top", legend = lgd, cex = 1.5, bty = "n")
    }
  }
  options(scipen = 0)
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
}

















sampleSize = nrow(datResv)




for (i in 1:length(sampleWeights))
{
  cat("Sample set index =", i, "\n")
  
  
  dat = datNonzeroResv
  ind = sample(nrow(dat), sampleSize, replace = F)
  dat = dat[ind, , drop = F]
  dat = dat[order(dat$MDR), , drop = F]
  
  
  system.time({p0andEmpDistrs = estimateEmpPMFs(
    mdr = dat$MDR, cdr = dat$CDR, windows, P0 = P0,
    targetMDRs = targetMDRs,
    regridMethod = "lmm",
    correctBias = FALSE,
    # correctBias = TRUE,
    nonDegenerateOldDistrs = NULL,
    # nonDegenerateOldDistrs = tohoku,
    # sampleWeightOnOldDistrs = 0.99,
    sampleWeightOnOldDistrs = sampleWeights[i],
    empDistrSupportSize = 64L,
    maxCore = 1000L,
    exponentialScalingForLowerMDRs = FALSE
  )})
  
  
  # Plot all empirical PMFs.
  if (F)
  {
    
    mdrInterest = unique(
      c(1:10, seq(10L, 100L, by = 10), seq(100L, 1000L, by = 100L),
        seq(1000L, 10000L, by = 1000L), seq(10000L, 100000L, by = 10000L))) / 1e5
    mdrInterest = mdrInterest[mdrInterest < 1]
    ind = match(as.integer(round(mdrInterest * 1e6)), 
                as.integer(round(p0andEmpDistrs$targetMDRs * 1e6)))
    tmp = p0andEmpDistrs$empDistrs[ind]
    ylim = c(0, max(unlist(lapply(tmp, function(x) 
      max(x[[2]][-c(1,length(x[[2]]))])))))
    
    
    pdf("../figure/tmpempDistrs.pdf", width = 13, height = 13 * (9/16))
    options(scipen = 999)
    par(mar = c(2, 3.3, 0, 0), mfrow = c(5, 9), family = "serif")
    for (i in 1:length(tmp))
    {
      if (i == 1) yaxt = 's' else yaxt = 'n'
      if (i == 37) xaxt = 's' else xaxt = 'n'
      h = wrapAsHist(tmp[[i]])
      plot(h, main = "",
           # type = "h", bty = "L", 
           xlim = c(0, 1), 
           ylim = ylim, cex.lab = 1.5, cex.axis = 1.5, xlab = "", col = "gray",
           ylab = "", las = 1, xaxt = xaxt, yaxt = yaxt, freq = F, border = NA)
      if (i == 1) title(xlab = "DR", line = 0.75, cex.lab = 1.5)
      if (i == 10) title(ylab = "P", line = 1, cex.lab = 1.5)
      if (i == 1) lgd = paste0("MDR = \n", mdrInterest[i])
      else lgd = paste0(mdrInterest[i])
      legend("top", legend = lgd, bty= "n", cex = 1.5)
    }
    options(scipen = 0)
    dev.off()
    
    
    pdf("../figure/tmpempDistrsMainIndividualScope.pdf", width = 13, height = 13 * (9/16))
    options(scipen = 999)
    par(mar = c(2, 3.3, 0, 0), mfrow = c(5, 9), family = "serif")
    for (i in 1:length(tmp))
    {
      if (i == 1) yaxt = 's' else yaxt = 'n'
      if (i == 37) xaxt = 's' else xaxt = 'n'
      h = wrapAsHist(tmp[[i]])
      ylim = c(0, max(h$density[-1]))
      plot(h, main = "",
           # type = "h", bty = "L", 
           # xlim = c(0, 1), 
           ylim = ylim,
           cex.lab = 1.5, cex.axis = 1.5, xlab = "", col = "gray",
           ylab = "", las = 1, xaxt = xaxt, yaxt = yaxt, freq = F, border = NA)
      if (i == 1) title(xlab = "DR", line = 0.75, cex.lab = 1.5)
      if (i == 10) title(ylab = "P", line = 1, cex.lab = 1.5)
      if (i == 1) lgd = paste0("MDR = \n", mdrInterest[i])
      else lgd = paste0(mdrInterest[i])
      legend("top", legend = lgd, bty= "n", cex = 1.5)
    }
    options(scipen = 0)
    dev.off()
    
    
    
  }
  
  
  # ==============================================================================
  # Create initial parameters.
  # ==============================================================================
  ndiv = 3L
  abcLB = c(1.01, 0.1, 0.1)
  abcUB = c(30, 30, 30)
  
  
  if (multipleInits)
  {
    abc = apply(rbind(abcLB, abcUB), 2, function(x)
    {
      exp(seq(log(x[1]), log(x[2]), len = ndiv + 2L)[2L:(2L + ndiv - 1L)])
    })
    abc = t(expand.grid(a = abc[,1], b = abc[,2], c = abc[,3]))
  }
  else 
  {
    abc = as.matrix(c(4, 5, 6))
  }
  
  
  tmp = extractMain(p0andEmpDistrs, targetMDRs, 
                    normalizeMainPart = TRUE, removeZeroPs = FALSE)
  condEmpDistrs = tmp$conditionalDistrs
  lm1 = tmp$lm1
  
  
  startingPmf = which.min(abs(unlist(lapply(
    p0andEmpDistrs, function(x) sum(x[[1]] * x[[2]]))) - startMDRtoFit))
  
  
  # Plot unnormalized main PMF
  if (F)
  {
    
    tmp = condEmpDistrs
    pdf("../figure/tmpempDistrsMain.pdf", width = 10, height = 10 * (9/16))
    par(mar = c(4.2, 5, 0, 0))
    for (i in 1:length(tmp))
    {
      if (i %% 1000L == 0L) cat(i, "")
      plot(as.data.frame(tmp[[i]]), type = "h", bty = "L", xlim = c(0, 1), 
           ylim = c(0, 0.1), cex.lab = 2, cex.axis = 1.5, xlab = "DR", 
           ylab = "P", las = 1)
      legend("top", legend = paste0("MDR = ", p0andEmpDistrs$targetMDRs[i]), 
             bty= "n", cex= 2)
    }
    dev.off()
  }
  
  
  # Optimization.
  if (T)
  {
    if (!sequentialUpdate) startingPmf = 1L
    
    
    secondHalf = startingPmf:length(condEmpDistrs)
    system.time({optRstSecondHalf = LBFGSBtrbFitList(
      abc = abc, 
      lm1 = lm1[secondHalf],
      empDistrList = condEmpDistrs[secondHalf], 
      abcLB = abcLB, abcUB = abcUB, 
      scaleEps = 1e-8, scaleMaxit = 100, 
      distanceFun = distanceFun,
      max_iterations = 100, maxCore = 1000, 
      RIBlib = "Numerical Recipes", 
      sequentialUpdate = sequentialUpdate, 
      hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
      epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
      max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
      ftol = 1e-4, wolfe = 0.9)
    })
    
    
    if (startingPmf != 1L)
    {
      firstHalf = startingPmf:1  
      system.time({optRstFirstHalf = LBFGSBtrbFitList(
        abc = abc, 
        lm1 = lm1[firstHalf],
        empDistrList = condEmpDistrs[firstHalf], 
        abcLB = abcLB, abcUB = abcUB, 
        scaleEps = 1e-8, scaleMaxit = 100, 
        distanceFun = distanceFun,
        max_iterations = 100, maxCore = 1000, 
        RIBlib = "Numerical Recipes", 
        sequentialUpdate = sequentialUpdate, 
        hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
        epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
        max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
        ftol = 1e-4, wolfe = 0.9)
      })
      optRst = list()
      optRst$param = cbind(optRstFirstHalf$param[, firstHalf], 
                           optRstSecondHalf$param[, -1])
      optRst$fval = c(optRstFirstHalf$fval[firstHalf], 
                      optRstSecondHalf$fval[-1])
      optRst$niter = c(optRstFirstHalf$niter[firstHalf], 
                       optRstSecondHalf$niter[-1])
    }
    else optRst = optRstSecondHalf
    
    
  }
  
  
  ensembleRst[[length(ensembleRst) + 1]] = list(
    optRst = optRst, empMeans = unlist(
      lapply(p0andEmpDistrs, function(x) keyALGs::Mean(x))))
}






























































plot(x = targetMDRs, y = P0, type = "l")
a = t(as.data.frame(lapply(ensembleRst, function(x)
{
  x$optRst$param[1,]
})))
dimnames(a) = NULL
a = apply(a, 2, function(x)
{
  mean(x[x >= abcLB[1] & x <= abcUB[1]])
})
# tmp = movingAverageSmoothing(a, 1001, 1)
plot(x = targetMDRs, y = a, type = "l", 
     ylim = c(abcLB[1], abcUB[1]),
     log = "x")



b = t(as.data.frame(lapply(ensembleRst, function(x)
{
  x$optRst$param[2,]
})))
dimnames(b) = NULL
b = apply(b, 2, function(x)
{
  mean(x[x >= abcLB[2] & x <= abcUB[2]])
})
plot(x = targetMDRs, y = b, type = "l", 
     # ylim = c(abcLB[2], abcUB[2])
     log = "x")


c = t(as.data.frame(lapply(ensembleRst, function(x)
{
  x$optRst$param[3,]
})))
dimnames(c) = NULL
c = apply(c, 2, function(x)
{
  mean(x[x >= abcLB[3] & x <= abcUB[3]])
})
plot(x = targetMDRs, y = c, type = "l", 
     # ylim = c(abcLB[3], abcUB[3])
     log = "x")


fval = t(as.data.frame(lapply(ensembleRstSeq, function(x) x$optRst$fval)))
fval = apply(fval, 2, function(x) mean(x[x > -100 & x < 100]))
plot(x = targetMDRs, y = fval, type = "l")






str(optRst)


d = solve_d(optRst$param, eps= 1e-8, maxit = 100)
abcd = rbind(optRst$param[1:3, , drop = F], d)
tailProbThreshold = 1e-4
maxes = pmin(1, actuar_qtrbeta(
  1 - tailProbThreshold / (1 - p0andEmpDistrs$P0), abcd))
ind = longestNonDecreasingSubseq(maxes)
maxFun = splinefun(x = p0andEmpDistrs$targetMDRs[ind], 
                   y = maxes[ind], method = "hyman")
maxes = maxFun(p0andEmpDistrs$targetMDRs)


inputParam = rbind(p0andEmpDistrs$targetMDRs, maxes, p0andEmpDistrs$P0, abcd)
distable = FTDFQA(
  inputParam, supportSize = 64, regridMethod = "lr", 
  RIBlib = "Numerical Recipes",
  outputProbsInRows = F, fineDiscretizationSize = 2000, maxCore = 1)


# Plot everything.
if (F)
{
  
  
  pdf("../figure/tmpAlldistrs.pdf", width = 10, height = 10 * (9/16))
  apply(distable, 2, function(x)
  {
    plot(x = seq(0, x[2], len = 64), y = x[-c(1,2)], xlim = c(0, 1), 
         ylim= c(0, 1), type = "h", bty = "L", xlab = "Damage ratio",
         ylab = "P", las = 1)
  })
  dev.off()
  
  
  pdf("../figure/tmpAlldistrsTransition.pdf", width = 10, height = 10 * (9/16))
  for (i in 3:nrow(distable))
  {
    if (i == 3) plot(x = distable[1,],  y = distable[i,], type = "l")
    else lines(x = distable[1,],  y = distable[i,], type = "l")
  }
  dev.off()
  
}














param = optRst$param
param[4, ] = targetMDRs
# ind = 1:ncol(param)
# ind = (ncol(param)-10L):ncol(param)
ind = 1:10
param = param[, ind, drop = F]
P0 = p0andEmpDistrs$P0[ind]
optRstDistrTable = makeDistrTable(
  param, 
  tentativeMax = 1, 
  lastProbLowerBound = 1e-6, 
  P0 = P0, 
  size = 64, 
  eps  = 1e-11, 
  maxit = 100, 
  maxCore = 1,
  useNewton = FALSE, 
  RIBlib = "Numerical Recipes",
  PfirstCoverFullDelta = T)


# Check problems.
if (F)
{
  
  Rcpp::sourceCpp("cpp/amalgamation.cpp")
  
  
  abc_lm1 = matrix(c(2, 1, 1, 1e-5))
  p0 = 0.98
  
  
  tmp = actuar::discretize(pgamma(x, 1), method = "unbiased",
                   lev = actuar::levgamma(x, 1), from = 0, to = 5, step = 0.5)
  
  
  rst = discretizeUnbiased(
    abc_lm1 = abc_lm1, 
    lastProbLowerBound = 1e-4,
    P0 = p0,
    size = 64, 
    eps = 1e-8, 
    maxit = 100,
    RIBlib = "R::pbeta")
  
  
  
  
  ind= 1L
  param = optRst$param[, ind]
  param = c(1, 1, 1, 1e-5)
  param[4] = targetMDRs[ind]
  P0 = p0andEmpDistrs$P0[ind]
  solve_d(param, eps = 1e-8, maxit = 100)
  qtrbeta(1-1e-4,param)
  # When b < 1, density at 0 is inifinity.
  
  
  rst = discretize(param, 
                   tentativeMax = 1,
                   lastProbLowerBound = 1e-3,
                   size = 1024^3, 
                   eps = 1e-8, 
                   maxit = 100,
                   useNewton = F, 
                   RIBlib = "Numerical Recipes",
                   PfirstCoverFullDelta = F)
  rst$distr$P
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  tmpParam = param[, 1]
  tmpP0 = P0[1]
  abc_lm1_ = c(tmpParam[1:3], tmpParam[4] / (1 - tmpP0))
  
  
  p  = 1e-3 / (1 - tmpP0)
  qtrbeta(1 - p, abc_lm1_ )
  
  
  distr = discretize(abc_lm1 = abc_lm1_, 
                     tentativeMax = 1,
                     lastProbLowerBound = p,
                     size = 63, 
                     eps = 1e-11, 
                     maxit = 100,
                     useNewton = F, RIBlib = "Numerical Recipes")
  
  
  distr = list(val = c(0, distr$distr$val), P = c(P0, distr$distr$P))
  distr$P[-1] = distr$P[-1] * (1 -P0)
  
  
  
}























system.time({optRst2 = LBFGSBtrbFitList(
  abc = abc, lm1 = lm1, empDistrList = condEmpDistrs, 
  abcLB = abcLB, abcUB = abcUB, 
  scaleEps = 1e-8, scaleMaxit = 100, distanceFun = "likelihood", 
  max_iterations = 100, maxCore = 1000, RIBlib = "Numerical Recipes", 
  sequentialUpdate = F, hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5, 
  epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10, 
  max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, 
  ftol = 1e-4, wolfe = 0.9)
})




rst = LBFGSBtrbFit(
  # c(abc[, 1], lm1[1]),
  c(abc[, 1], lm1[1]),
  # c(2,3, 4,0.5),
  empDistr = empDistrs[[1]], 
  abcLB = c(1.001, 0.001, 0.1), 
  abcUB = c(30, 30, 30),
  scaleEps = 1e-8, scaleMaxit = 100, distanceFun = "likelihood",
  max_iterations = 100, hgrad = 0, centralDiff = T, m = 6, epsilon = 1e-5,
  epsilon_rel = 1e-5, past = 1, delta = 1e-10, max_submin = 10,
  max_linesearch = 20, min_step = 1e-20, max_step = 1e+20, ftol = 1e-4,
  wolfe = 0.9); rst; abc_lm1; init


a = log(2); b = log(3); c = log(4); lm1 = log(0.5)
s = 0.01
N = 1000L
lb = c(log(1.001), log(0.01), log(1), log(1e-10))
ub = c(log(30), log(30), log(30), log(1))
params = list(c(a, b, c, lm1))
for ( i in 2:N)
{
  tmp = params[[i - 1]] + rnorm(4, 0, s)
  params[[i]] = pmax(lb, pmin(tmp, ub))
}
abc = exp(matrix(unlist(params), nrow = 4))
lm1 = abc[4, , drop = T]
abc = abc[1:3, , drop = F]
abcLB = exp(lb)
abcUB = exp(ub)



LBFGSBtrbFitList(
  NumericMatrix abc,
  NumericVector lm1,
  List empDistrList,
  NumericVector abcLB = NumericVector(0),
  NumericVector abcUB = NumericVector(0),
  double scaleEps = 1e-8, 
  int scaleMaxit = 100,
  String distanceFun = "likelihood",
  int max_iterations = 100,
  int maxCore = 15,
  String RIBlib = "Numerical Recipes",
  bool sequentialUpdate = false,
  double hgrad = 0,
  bool centralDiff = true,
  int m = 6,
  double epsilon = 1e-5,
  double epsilon_rel = 1e-5,
  int past = 1,
  double delta = 1e-10,
  int max_submin = 10,
  int max_linesearch = 20,
  double min_step = 1e-20,
  double max_step = 1e+20,
  double ftol = 1e-4,
  double wolfe = 0.9
)






