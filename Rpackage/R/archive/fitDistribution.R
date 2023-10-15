
# ==============================================================================
# Read data in.
# ==============================================================================
load("../data/cvgAdata.Rdata")
dat = cvgData[order(cvgData$MDR, cvgData$CDR), ]
rm(cvgData); gc()
# ==============================================================================


# A claim’s damage ratio is considered to be zero if it is <= 0.002
dat$CDR[dat$CDR <= 0.002] = 0
source('R/rfuns.R')


windowSize = 1000L
windowSpeed = 30L
wdws = computeWindows(nrow(dat), windowSize = windowSize, 
                      speed = windowSpeed)


mvAvgs = lapply(as.integer(c(125, 300, 1000, 3000)), function(x)
{
  wds = computeWindows(nrow(dat), windowSize = x, 
                       speed = windowSpeed)
  
  
  avgMdr = unlist(lapply(wds, function(x)
  {
    mean(dat$MDR[x[1]:x[2]])
  }))
  
  
  avgCdr = unlist(lapply(wds, function(x)
  {
    mean(dat$CDR[x[1]:x[2]])
  }))
  
  
  avgP0 = unlist(lapply(wds, function(x)
  {
    sum(dat$CDR[x[1]:x[2]] <= 0) / (x[2] - x[1] + 1)
  }))
  
  
  data.frame(mdr = avgMdr, cdr = avgCdr, p0 = avgP0)
})


names(mvAvgs) = c('wsize125', 'wsize300', 'wsize1000', 'wsize3000')


# Print some stats.
if (T)
{
  
  mdrRange = range(unlist(lapply(mvAvgs, function(x) range(x$mdr))))
  cdrRange = c(0, max(unlist(lapply(mvAvgs, function(x) max(x$cdr)))))
  p0Range = c(0, max(unlist(lapply(mvAvgs, function(x) max(x$p0)))))
  dir.create('../figure', showWarnings = F)
  pdf('../figure/mvAvgsDifferengWindowSizes.pdf', width = 8, height = 8)
  par(mar = c(4, 5, 0, 0), family = "serif", mfrow = c(4, 2))
  wsize = c(125, 300, 1000, 3000)
  for (i in 1:length(mvAvgs))
  {
    if (i == 1) ylab = "CDR" else ylab = ""
    if (i == 4) xlab = "MDR" else xlab = ""
    plot(x = mvAvgs[[i]]$mdr, y = mvAvgs[[i]]$cdr, type = "l", 
         col = scales::alpha('red', 0.5), 
         cex.lab = 2, cex.axis = 1.5, bty = "L", ylab = ylab, 
         xlim = mdrRange, ylim = cdrRange, xlab = xlab, las = 1)
    
    
    if (i == 1) legend(
      "bottomright", legend = c(expression(CDR), expression(P[0])),
      cex = 2, bty = "n", col = scales::alpha(c('red', 'blue'), 0.5), 
      lwd = c(2, 2))
    
    
    if (i == 1) ylab = expression(P[0]) else ylab = ""
    plot(x = mvAvgs[[i]]$mdr, y = mvAvgs[[i]]$p0, type = "l", 
         col = scales::alpha('blue', 0.5), 
         cex.lab = 2, cex.axis = 1.5, bty = "L", ylab = ylab, 
         xlim = mdrRange, ylim = p0Range, xlab = "", las = 1)
    
    
    legend(
      "topright", 
      legend = c(paste0("Window size = ", wsize[i]), 
                 paste0('Sliding speed = ', 30)), 
      bty = "n",
      cex = 2)
  }
  dev.off()
  
  
}


# Window size = 1% * data size, sliding speed = 5% * 1% * data size
windowSize = as.integer(round(0.01 * nrow(dat))) # 1784
windowSpeed = as.integer(round(0.01 * 0.01 * nrow(dat))) # 18
wdws = computeWindows(nrow(dat), windowSize = windowSize, 
                      speed = windowSpeed)
mvAvgs = as.data.frame(t(matrix(unlist(lapply(wdws, function(x)
{
  c(mean(dat$MDR[x[1]:x[2]]), mean(dat$CDR[x[1]:x[2]]), sum(
    dat$CDR[x[1]:x[2]] <= 0) / (x[2] - x[1] + 1L))
})), nrow = 3)))
colnames(mvAvgs) = c("MDR", "CDR", "P0")


# Plot name: window size 1% of data size; sliding speed 5% of window size.
pdf('../figure/ws-1pct-speed-1pct.pdf', width = 10, height = 10 * 0.618) 
par(mar = c(4, 5, 0, 0), family = 'serif')
plot(x = mvAvgs$MDR, y = mvAvgs$P0, bty = "L", 
     col = 'darkblue', type = "l", xlab = "MDR", ylab = expression(P[0]),
     cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(0, 1), ylim = c(0, 1))
lines(x = c(0, mvAvgs$MDR[1]), y = c(1, mvAvgs$P0[1]), 
      col = "red")
lines(x = c(mvAvgs$MDR[nrow(mvAvgs)], 1), 
      y = c(mvAvgs$P0[nrow(mvAvgs)], 0), 
      col = "red")
legend("topright", legend = c("Data", "Linked"), lwd = c(1, 1), bty = "n",
       cex = 2, col = c("darkblue", "red"))
legend("right", legend = c(
  expression('Window size = 1784 (1% of data size)'), 
  expression('Sliding speed = 18 (1% of window size)'),
  expression('N(window)'%~~%10000)
  ), 
  bty = "n", cex = 2)
legend("top", legend = expression("MDR vs."~P[0]~"moving average"), 
       cex = 2, bty = "n")
dev.off()


# ==============================================================================
# Longest decreasing subsequence.
# ==============================================================================
source('R/sourceCppCharlie.R')
sourceCppCharlie('cpp/longestSubseq.cpp', verbose = T)
mdr = c(0, mvAvgs$MDR, 1)
p0 = c(1, mvAvgs$P0, 0)
lds = longestIncreasingSubseq(-p0)
p0 = p0[lds]
mdr = mdr[lds]


# ==============================================================================
# Plot data against the longest subsequence.
# ==============================================================================
pdf('../figure/ws-1pct-speed-1pct-lds.pdf', width = 10, height = 10 * 0.618) 
par(mar = c(4, 5, 0, 0), family = 'serif')
plot(x = mvAvgs$MDR, y = mvAvgs$P0, bty = "L", 
     col = 'darkblue', type = "l", xlab = "MDR", ylab = expression(P[0]),
     cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(0, 1), ylim = c(0, 1))
lines(x = c(0, mvAvgs$MDR[1]), y = c(1, mvAvgs$P0[1]), col = "red")
lines(x = c(mvAvgs$MDR[nrow(mvAvgs)], 1), 
      y = c(mvAvgs$P0[nrow(mvAvgs)], 0), 
      col = "red")
lines(x = mdr, y = p0, type = "p", col = "olivedrab3", pch = 16, cex = 0.75)
legend("right", legend = c("Data", "Linked", "Longest decreasing\nsubsequence"), 
       lwd = c(1, 1, NA), bty = "n",
       cex = 2, col = c("darkblue", "red", "olivedrab3"), pch = c(NA, NA, 16))
legend("top", legend = expression("MDR vs."~P[0]~"moving average"), 
       cex = 2, bty = "n")
dev.off()


# ==============================================================================
# Plot linear interpolated P0
# ==============================================================================
lrfun = approxfun(x = mdr, y = p0)
sampledMDR = c(mdr[1], seq(mdr[2], mdr[length(mdr)], len = 300))
sampledP0 = lrfun(sampledMDR)
pdf('../figure/ws-1pct-speed-1pct-lds-linear.pdf', width = 10, height = 10 * 0.618) 
par(mar = c(4, 5, 0, 0), family = 'serif')
plot(x = mvAvgs$MDR, y = mvAvgs$P0, bty = "L", 
     col = 'darkblue', type = "l", xlab = "MDR", ylab = expression(P[0]),
     cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(0, 1), ylim = c(0, 1))
lines(x = c(0, mvAvgs$MDR[1]), y = c(1, mvAvgs$P0[1]), col = "red")
lines(x = c(mvAvgs$MDR[nrow(mvAvgs)], 1), 
      y = c(mvAvgs$P0[nrow(mvAvgs)], 0), 
      col = "red")
lines(x = sampledMDR, y = sampledP0, type = "l", col = "olivedrab3", lwd = 2)
legend("right", legend = c("Data", "Linked", "Linear interpolation"), 
       lwd = c(1, 1, 2), bty = "n",
       cex = 2, col = c("darkblue", "red", "olivedrab3"))
legend("top", legend = expression("MDR vs."~P[0]~"moving average"), 
       cex = 2, bty = "n")
dev.off()


# ==============================================================================
# Plot Hyman interpolated P0
# ==============================================================================
cbfun = splinefun(x = mdr, y = p0, method = 'hyman')
sampledMDR = c(mdr[1], seq(mdr[2], mdr[length(mdr)], len = 300))
sampledP0 = cbfun(sampledMDR)
pdf('../figure/ws-1pct-speed-1pct-lds-hyman.pdf', width = 10, height = 10 * 0.618) 
par(mar = c(4, 5, 0, 0), family = 'serif')
plot(x = mvAvgs$MDR, y = mvAvgs$P0, bty = "L", 
     col = 'darkblue', type = "l", xlab = "MDR", ylab = expression(P[0]),
     cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(0, 1), ylim = c(0, 1))
lines(x = c(0, mvAvgs$MDR[1]), y = c(1, mvAvgs$P0[1]), col = "red")
lines(x = c(mvAvgs$MDR[nrow(mvAvgs)], 1), 
      y = c(mvAvgs$P0[nrow(mvAvgs)], 0), 
      col = "red")
lines(x = sampledMDR, y = sampledP0, type = "l", col = "olivedrab3", lwd = 2)
legend("right", legend = c("Data", "Linked", "Hyman interpolation"), 
       lwd = c(1, 1, 2), bty = "n",
       cex = 2, col = c("darkblue", "red", "olivedrab3"))
legend("top", legend = expression("MDR vs."~P[0]~"moving average"), 
       cex = 2, bty = "n")
dev.off()


# ==============================================================================
# Take samples from interpolated function, etc.
# ==============================================================================
sampledMDR = c(mdr[1], seq(mdr[2], mdr[length(mdr)], len = 300))
sampledP0 = lrfun(sampledMDR)
lrfunInv = approxfun(x = p0, y = mdr)
sampledMDRfromInterpolation = c(seq(0, 1, len = 100), 
               lrfunInv(seq(0, 1, len = 100)))
sampledMDRfromInterpolation = sort(unique(sampledMDRfromInterpolation))
sampledP0fromInterpolation = lrfun(sampledMDRfromInterpolation)
pdf('../figure/ws-1pct-speed-1pct-lds-linear-takeSamples.pdf', 
    width = 10, height = 10 * 0.618) 
par(mar = c(4, 5, 0, 0), family = 'serif')
plot(x = mvAvgs$MDR, y = mvAvgs$P0, bty = "L", 
     col = 'darkblue', type = "l", xlab = "MDR", ylab = expression(P[0]),
     cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(0, 1), ylim = c(0, 1))
lines(x = c(0, mvAvgs$MDR[1]), y = c(1, mvAvgs$P0[1]), col = "red")
lines(x = c(mvAvgs$MDR[nrow(mvAvgs)], 1), 
      y = c(mvAvgs$P0[nrow(mvAvgs)], 0), 
      col = "red")
lines(x = sampledMDR, y = sampledP0, type = "l", col = "olivedrab3", lwd = 2)
lines(x = sampledMDRfromInterpolation, y = sampledP0fromInterpolation,
      col = "black", type = "p")
legend("right", legend = c("Data", "Linked", "Linear interpolation", " ", 
                           "Uniform samples from\ninterpolated function"), 
       lwd = c(1, 1, 2, NA, NA), bty = "n",
       cex = 2, col = c("darkblue", "red", "olivedrab3", "white", "black"), 
       pch = c(NA, NA, NA, NA, 1))
legend("top", legend = expression("MDR vs."~P[0]~"moving average"), 
       cex = 2, bty = "n")
dev.off()




# ==============================================================================
# Fit a function in the form of transformed beta survival function.
# ==============================================================================
# exp(-b) * (1 - x ^ a * exp(b * x ^ c)) + (1 - exp(-cfs[2]))
fmdrp0 = function(x, cfs)
{
  exp(-cfs[2]) * (1 - x ^ cfs[1] * exp(cfs[2] * x ^ cfs[3])) + 
    (1 - exp(-cfs[2]))
}


# Tuning code.
if (F)
{
  plot(x = sampledMDRfromInterpolation, y = sampledP0fromInterpolation)
  lines(x = sampledMDRfromInterpolation, 
        y = fmdrp0(sampledMDRfromInterpolation, 
                   c(0.05, 0.08, 1)), col = "red")  
}


# ==============================================================================
# Take samples from interpolated function, etc.
# ==============================================================================
sampledMDR = c(mdr[1], seq(mdr[2], mdr[length(mdr)], len = 300))
sampledP0 = lrfun(sampledMDR)
lrfunInv = approxfun(x = p0, y = mdr)
sampledMDRfromInterpolation = c(seq(0, 1, len = 100), 
                                lrfunInv(seq(0, 1, len = 100)))
sampledMDRfromInterpolation = sort(unique(sampledMDRfromInterpolation))
sampledP0fromInterpolation = lrfun(sampledMDRfromInterpolation)
pdf('../figure/ws-1pct-speed-1pct-lds-linear-monoConvexFit.pdf', 
    width = 10, height = 10 * 0.618) 
par(mar = c(4, 5, 0, 0), family = 'serif')
plot(x = mvAvgs$MDR, y = mvAvgs$P0, bty = "L", 
     col = 'darkblue', type = "l", xlab = "MDR", ylab = expression(P[0]),
     cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(0, 1), ylim = c(0, 1))
lines(x = c(0, mvAvgs$MDR[1]), y = c(1, mvAvgs$P0[1]), col = "red")
lines(x = c(mvAvgs$MDR[nrow(mvAvgs)], 1), 
      y = c(mvAvgs$P0[nrow(mvAvgs)], 0), 
      col = "red")
lines(x = sampledMDR, y = sampledP0, type = "l", col = "olivedrab3", lwd = 2)
lines(x = sampledMDRfromInterpolation, y = sampledP0fromInterpolation,
      col = "black", type = "p")
lines(x = sampledMDR, y = fmdrp0(sampledMDR, c(0.05, 0.08, 1)), 
      col = scales::alpha('orange', 0.7), lwd = 4)
legend("right", legend = c("Data", "Linked", "Linear interpolation", 
                           "Uniform samples from\ninterpolated function",
                           "A monotone convex fit"), 
       lwd = c(1, 1, 2, NA, 3), bty = "n",
       cex = 2, col = c("darkblue", "red", "olivedrab3", "black", "orange"), 
       pch = c(NA, NA, NA, 1, NA))
legend("top", legend = expression("MDR vs."~P[0]~"moving average"), 
       cex = 2, bty = "n")
dev.off()



allmdrs = unique(as.integer(round(c(
  seq(0, 0.1, by = 1e-5), seq(0.1, 1, by = 1e-4)) * 1e6))) / 1e6
allmdrs = allmdrs[-c(1, length(allmdrs))]
allp0 = lrfun(allmdrs)
allCondMdrs = allmdrs / (1 - allp0)
mdrP0condMdr = data.frame(mdr = allmdrs, p0 = allp0, condMdr = allCondMdrs)
rm(allmdrs, allp0, allCondMdrs); gc()


findCondMDR = function(mdr, thedata = mdrP0condMdr)
{
  inds = match(as.integer(round(mdr * 1e5)), 
               as.integer(round(thedata$mdr * 1e5)))
  thedata$allCondMdrs[inds]
}


# Test MDR = 0 05
roundedMDRs = unlist(lapply( wdws, function(x)
{
  m = mean(dat$MDR[x[1]:x[2]])
  if (m < 0.1) round(m, 5)
  else round(m, 4)
}))
tmp = aggregate(list(1:length(wdws)), list(roundedMDRs = roundedMDRs), 
                function(x) x, simplify = F)
tmp[[2]] = lapply(tmp[[2]], function(x)
{
  range(unlist(wdws[x]))
})
mdrWdws = as.list(tmp)
names(mdrWdws) = c("mdr", "windowBounds")
rm(tmp); gc()


tmp = match(mdrWdws$mdr, mdrP0condMdr$mdr)
mdrWdws$p0 = mdrP0condMdr$p0[tmp]
mdrWdws$condMdr = mdrP0condMdr$condMdr[tmp]


# tmp = mdrWdws
mdrWdws$cdr = dat$CDR
mdrWdws$windowBounds = matrix(unlist(mdrWdws$windowBounds), nrow = 2)
mdrWdws$retriveCDRandTargetMean = function(m)
{
  k = which(mdrWdws$mdr == m)
  b = mdrWdws$windowBounds[, k]
  r = mdrWdws$cdr[b[1]:b[2]]
  list(mdr = m, condMdr = mdrWdws$condMdr[k], nonzeroCdrs = r[r > 0],
       p0 = mdrWdws$p0[k])
}


testDat = mdrWdws$retriveCDRandTargetMean(0.05)












































mdrWdws$p0 = tmp$p0
mdrWdws$condMdr = tmp$condMdr
mdrWdws$windowBounds = matrix(unlist(tmp$windowBounds), nrow = 2)
mdrWdws$mdr = tmp$mdr
mdrWdws$retriveCDRandTargetMean = function(m)
{
  
}






























