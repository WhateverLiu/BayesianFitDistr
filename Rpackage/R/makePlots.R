

pmfTableToDistrList = function(X)
{
  npoint = nrow(X) - 2L
  apply(X, 2, function(x)
  {
    list(val = seq(0, x[2], len = npoint), P = x[-c(1, 2)])
  })
}


#'
#' Plot \eqn{P_0}
#' 
#' Plot fitted \eqn{P_0} versus MDRs.
#' 
#' @param figPath  Path to the figure to be saved. File extension MUST be
#' \code{.png} or \code{.PNG}.
#' 
#' @param targetMDRs  A numeric vector of MDRs.
#' 
#' @param p0models  A list returned from function \code{\link{ensembleP0}}.
#' 
#' @return None
#' 
#' @details  A helper function for visualization. Please walk 
#' through \code{vignette('Workflow', package = 'NGFMfitDistr')} to see how it 
#' is used.
#'
plotP0model = function(figPath, targetMDRs, p0models)
{
  P0 = rowMeans(as.data.frame(p0models))
  png(figPath, width = 10, height = 0.618 * 10, units = "in", res = 120)
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




#'
#' Plot TrB Parameters
#' 
#' Plot trained Transformed Beta parameters versus MDRs.
#' 
#' @inheritParams  plotP0model
#' 
#' @param optRst  Optimization result from function \code{\link{LBFGSBtrbFitList}}.
#' 
#' @param startingPmf  An integer as the index of the first PMF that was fitted
#' during bi-directional sequential fitting.
#' 
#' @param abc  A numeric vector of size 3. The initial \code{a, b, c} for fitting
#' the \code{startingPmf}-th empirical PMF.
#' 
#' @return None
#' 
#' @inherit plotP0model details
#' 
#'
plotTrBparameters = function(
    figPath, 
    optRst,
    startingPmf, 
    abc
    )
{

  
  png(figPath, width = 9, height = 9 * 0.618, units = "in", res = 120)
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
  text(x = xcoor, y = ycoor, labels = paste0(
    "Initial a b c = ", paste0(abc, collapse = " ")), cex = 2)
  dev.off()
  
}




wrapAsHist = function(X)
{
  val = X$val
  m = (val[-1] + val[-length(val)]) / 2
  sp = c(val[1] - (m[1] - val[1]), m, 
         val[length(val)] + (val[length(val)] - m[length(m)]))
  H = hist(val, breaks = sp, plot = F)
  H$density = X$P; H
}




#'
#' Plot 45 PMFs in PMF table
#' 
#' Plot PMFs associated with MDR = 0.00001, 0.00002, ..., 0.0001, 0.0002, ...
#' 0.001, 0.002, ..., 0.01, 0.02, ..., 0.1, 0.2, ..., 0.9.
#' 
#' @param figPath  Path to the figure to be saved. File extension MUST be
#' \code{.pdf}.
#' 
#' @param pmfTable  The distribution table as a matrix. Every column is 
#' \code{MDR, max, P0, P1, ..., Pmax}.
#' 
#' @param empDistrsLists  Result from \code{\link{estimateEmpPMFs}}.
#' 
#' @param mdrRangeInData  A numeric vector of size 2. Range of MDRs whose 
#' associated empirical PMFs are inferred from data. 
#' 
#' @return None
#' 
#' @inherit plotP0model details
#'
#'
plot45pmfsInTable = function(
    figPath,
    pmfTable,
    empDistrsLists,
    mdrRangeInData
  )
{
  tmp = unique(as.integer(round(c(
    seq(1e-5, 1e-4, by = 1e-5), 
    seq(1e-4, 1e-3, by = 1e-4), 
    seq(1e-3, 1e-2, by = 1e-3), 
    seq(1e-2, 1e-1, by = 1e-2), 
    seq(1e-1, 1-1e-1, by = 1e-1)) * 1e5)))
  
  
  targetMDRs = pmfTable[1, , drop = TRUE]
  whichEmpDistrs = match(tmp, as.integer(round(targetMDRs * 1e5)))
  
  
  empDistrsToPlot = empDistrsLists$biased[whichEmpDistrs]
  empDistrsBiasCorrectedToPlot = empDistrsLists$biasCorrected[whichEmpDistrs]
  
  
  whichPMFinTable2plot = match(tmp, as.integer(round(pmfTable[1,] * 1e5)))
  tablePmf2plot = pmfTableToDistrList(
    pmfTable[, whichPMFinTable2plot, drop = F])
  
  
  pdf(figPath, width = 10, height = 10 * 0.618)
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
    if (Xhist$mids[l] >= 1 - 1e-10) 
      Xhist$density[-c(1, l)] = Xhist$density[-c(1, l)] * s
    else Xhist$density[-c(1)] = Xhist$density[-c(1)] * s
    xlim = c(0, max(c(range(Y$val), range(Xhist$breaks), range(Z$val))))
    ylim = range(c(0, max(Y$P), max(Xhist$density), max(Z$P)))
    
    
    hcol = "gray"
    plot(Xhist, freq = FALSE, border = NA, col = hcol, main = "", 
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
      tcol = c("black", "black")
    }
    else  {
      lgd = c("Fitted & discretized", 
              "Extrapolated empirical. Fake data as surrogate!")
      tcol = c("black", "red")
    }
    legend("center", legend = lgd, bty = "n", cex = 1.5, 
           col = c("black", "gray"), pch = c(NA, 15), lwd = c(1.5, NA), 
           text.col = tcol)
    
    
    Zhist = wrapAsHist(Z)
    s = (diff(range(Y$val)) / (length(Y$val) - 1)  ) / 
      (diff(range(Z$val)) / (length(Z$val) - 1) )
    l = length(Zhist$mids)
    if (Zhist$mids[l] >= 1 - 1e-10) 
      Zhist$density[-c(1, l)] = Zhist$density[-c(1, l)] * s
    else Zhist$density[-c(1)] = Zhist$density[-c(1)] * s
    plot(Zhist, freq = FALSE, border = NA, col = hcol, main = "", 
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




#'
#' Plot 207 PMFs in PMF table
#' 
#' Plot PMFs associated with MDR = 0.00001, 0.00002, ..., 0.0001, 0.0002, ...
#' 0.001, 0.002, ..., 0.01, 0.02, ..., 0.1, 0.2, ..., 0.9.
#' 
#' @param figPath  Path to the figure to be saved. File extension MUST be
#' \code{.pdf}.
#' 
#' @param pmfTable  The distribution table as a matrix. Every column is 
#' \code{MDR, max, P0, P1, ..., Pmax}.
#' 
#' @param empDistrsLists  Result from \code{\link{estimateEmpPMFs}}.
#' 
#' @inheritParams plot45pmfsInTable
#' 
#' @return None
#' 
#' @inherit plotP0model details
#'
#'
plot207pmfsInTable = function(
    figPath,
    pmfTable,
    empDistrsLists,
    mdrRangeInData
)
{
  tmp = unique(as.integer(round(c(
    seq(1e-5, 1e-4, by = 1e-5), seq(1e-4, 1e-2, by = 1e-4), 
    seq(1e-2, 1 - 1e-2, by = 1e-2)) * 1e5)))
  targetMDRs = pmfTable[1, , drop = TRUE]
  whichEmpDistrs = match(tmp, as.integer(round(targetMDRs * 1e5)))
  empDistrsToPlot = empDistrsLists$biased[whichEmpDistrs]
  empDistrsBiasCorrectedToPlot = empDistrsLists$biasCorrected[whichEmpDistrs]
  
  
  whichPMFinTable2plot = match(tmp, as.integer(round(pmfTable[1,] * 1e5)))
  tablePmf2plot = pmfTableToDistrList(pmfTable[, whichPMFinTable2plot, drop = F])
  
  
  pdf(figPath, width = 10, height = 10 * 0.618)
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
    plot(Xhist, freq = FALSE, border = NA, col = hcol, main = "", 
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
      tcol = c("black", "black")
    }
    else  {
      lgd = c("Fitted & discretized", 
              "Extrapolated empirical. Fake data as surrogate!")
      tcol = c("black", "red")
    }
    legend("center", legend = lgd, bty= "n", cex = 1.5, 
           col = c("black", "gray"), pch = c(NA, 15),
           lwd = c(1.5, NA), text.col = tcol)
    
    
    Zhist = wrapAsHist(Z)
    s = (diff(range(Y$val)) / (length(Y$val) - 1)  ) / 
      (diff(range(Z$val)) / (length(Z$val) - 1) )
    l = length(Zhist$mids)
    if (Zhist$mids[l] >= 1-1e-10) 
      Zhist$density[-c(1, l)] = Zhist$density[-c(1, l)] * s
    else Zhist$density[-c(1)] = Zhist$density[-c(1)] * s
    plot(Zhist, freq = FALSE, border = NA, col = hcol, main = "", 
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
  



#'
#' Visualize Parameter Trends in PMF table
#' 
#' Plot parameter trends in QA's interest, namely \eqn{P_0}s, \eqn{P_{\max}}s, 
#' maxes, PMF means, coefficients of variation against MDR in a PMF table.
#' 
#' @param pmfTable  PMF table as a numeric matrix. Every
#' column is a (MDR, max, P0, P1, ..., Pmax). The first and last PMFs will be 
#' excluded from visualization.
#' 
#' @param figPath  Path to the figure to be saved. File extension must be
#' \code{.png} or \code{.PNG}.
#' 
#' @return None
#' 
plotQAparams = function(pmfTable, figPath)
{
  png(figPath, width = 10, height = 10 * 0.618, units = "in", res = 120)
  par(mar = c(4.2, 5, 0, 0), family = "serif")
  layout(rbind(c(1, 1, 2, 2, 3, 3),
               c(4, 5, 5, 6, 6, 7)))
  exclu = -c(1L, ncol(pmfTable))
  xs = pmfTable[1, exclu, drop = FALSE]
  plot(x = xs, y = pmfTable[2, exclu], type = "l", 
       col = "darkblue", bty = "L", xlim = range(xs), ylim = c(0, 1),
       ylab = "Max", cex.lab = 2, cex.axis = 1.5, las = 1, xlab = "" )
  plot(x = xs, y = pmfTable[3, exclu], col = "darkblue", type = "l", 
       ylab = expression(P[0]), cex.lab = 2, cex.axis = 1.5, bty = "L", 
       las = 1, xlab = "" )
  plot(x = xs, y = pmfTable[nrow(pmfTable), exclu], 
       col = "darkblue", cex.lab = 2, cex.axis = 1.5, bty = "L", 
       las = 1, ylab = expression(P[max]), type = "l", xlab = "" )
  plot.new()
  ms = apply(pmfTable[, exclu], 2, function(x)
  {
    sum(seq(0, x[2], len = nrow(pmfTable) - 2) * x[-c(1, 2)])
  })
  plot(x = xs, y = ms, col = "darkblue", cex.lab = 2, cex.axis = 1.5, 
       bty = "L", las = 1, ylab = expression(mu), type = "l", xlab = "MDR" )
  sds = sqrt(apply(pmfTable[, exclu], 2, function(x)
  {
    sum(seq(0, x[2], len = nrow(pmfTable) - 2) ^ 2 * x[-c(1, 2)])
  }) - ms ^ 2)
  cv = sds / ms
  plot(x = xs, y = cv, col = "darkblue", cex.lab = 2, cex.axis = 1.5, 
       bty = "L", las = 1, ylab = expression(sigma/mu), type = "l", xlab = "" )
  dev.off()  
}

