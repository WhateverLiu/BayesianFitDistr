

# scoringTable = pmftablelist$tuned$p64
# save(scoringTable, file = "../data/scoringTable.RData")
scoringTable = data.table::setDF(  data.table::fread(
  "../tempFiles/CovA_nPts64-debora.csv"))
scoringTable = t(scoringTable)
dimnames(scoringTable) = NULL
scoringTable = scoringTable[, -c(1, ncol(scoringTable))]


sp = seq(0, 1, len = 64)
tmp = apply(scoringTable, 2, function(x)
{
  if (x[2] >= 1) x[-c(1, 2)]
  else
  {
    X = list(val = seq(0, x[2], len = 64), P = x[-c(1, 2)])
    Y = keyALGs::rglr(X, ngd = sp) 
    Y[[2]]
  }
})


X = tmp
X = t(tmp)
Y = rank(X)
Y = matrix((Y - min(Y)) / diff(range(Y)), nrow = nrow(X))
quat = 5e-3
k = (1 - quat - quat) / diff(range(Y))
b = quat - k * min(Y)
Y = k * Y + b
Y = qnorm(Y)
Y = Y - min(Y)
Y = Y[c(seq(10L, 10000L, by = 10L), 10001:nrow(Y)) , ]


clr = unique(colorRampPalette(c("blue",   "white",  "red"))(10000))
clrSmall = unique(colorRampPalette(c("blue",   "white",  "red"))(100))


png("../figure/distTableImage.png", res = 120, width = 7, height = 7 / 1.3, 
    units = "in")
par(mar = c(4.2, 5, 0, 0), family = "serif")
image( 
  t(Y[nrow(Y):1, ]) [, nrow(Y):1], 
  bty = "n", xlab = "Damage ratio", ylab = "MDR",
  cex.lab = 2, cex.axis = 1.5, col = clr, las = 1, useRaster = T, 
  xlim = c(0, 1.3), xaxt = "n")
axis(side = 1, at = seq(0, 1, len = 5), labels = seq(0, 1, len = 5), 
     cex.axis = 1.5)
lines(x = rep(1.15, 100), y = seq(0.3, 0.7, len = 100), pch = 15, 
      col = clrSmall, type = "p")
text(x = 1.175, y = 0.275, labels = "Low Prob", cex = 1.25)
text(x = 1.175, y = 0.725, labels = "High Prob", cex = 1.25)
dev.off()


png("../figure/distTableImageWithClaims.png", res = 120, width = 7, height = 7 / 1.3, 
    units = "in")
par(mar = c(4.2, 5, 0, 0), family = "serif")
image( 
  t(Y[nrow(Y):1, ]) [, nrow(Y):1],
  bty = "n", xlab = "Damage ratio", ylab = "MDR",
  cex.lab = 2, cex.axis = 1.5, col = clr, las = 1, useRaster = T, xlim = c(0, 1.3),
  xaxt = "n")
axis(side = 1, at = seq(0, 1, len = 5), labels = seq(0, 1, len = 5), cex.axis = 1.5)
lines(x = rep(1.15, 100), y = seq(0.3, 0.7, len = 100), pch = 15, 
      col = clrSmall, type = "p")
text(x = 1.175, y = 0.275, labels = "Low Prob", cex = 1.25)
text(x = 1.175, y = 0.725, labels = "High Prob", cex = 1.25)
tmpClaims = datResv[c("MDR", "CDR")]
tmpClaims = unique(tmpClaims)
lines(tmpClaims, cex = 0.1, type= "p", col = scales::alpha("black", 0.5))
legend("topright", legend = "Claims", bty = "n", cex = 1.5, pch = 16)
dev.off()





tmpsp = seq(0, by = 0.05, len = 10)
tmp = dgamma(tmpsp, shape = 2, scale = 0.2)
tmp = tmp /sum(tmp)
pdf("../figure/pointWithinPmfScore.pdf", width = 10, height = 10 * 0.618)
par(mar = c(4.2, 5, 0, 0), family = "serif")
plot(x = tmpsp, y = tmp, type = "h", xlab = "Damage ratio", ylab = "P", 
     bty = "L", cex.lab = 2, cex.axis = 1.5, las = 1, ylim = c(0, max(tmp) * 1.25),
     xlim = c(0, 0.7))
axis(side = 1, at = tmpsp[length(tmpsp)], 
     labels = latex2exp::TeX('$x_{\\max}$'), cex.axis = 1.5)
text(x = tmpsp[length(tmpsp)], y = tmp[length(tmp)] * 1.1, labels = 
       latex2exp::TeX('$p_{\\max}$'), cex = 1.5)
text(x = 0.075, y = 0.001, labels = latex2exp::TeX('$\\Delta x$'), cex = 1.5 )
mid = mean(tmpsp[7:8])
pmid = mean(tmp[7:8])
lines(x = mid, y = 0, type = "p", cex = 1.5, pch = 16, col = "blue")
lines(x = mid, y = pmid, type = "h", lty = 2, col = "blue")
lines(x = tmpsp[7:8], y = tmp[7:8], lty = 2, col = "blue")
lines(x = 0.6, y = 0, type = "p", cex = 1.5, pch = 16, col = "red")
text(x = 0.6, y = 0.008, labels = "?", cex = 1.5)
legend("topleft", legend = c("u: claim within PMF support", 
                             "v: claim out of PMF support"),
       pch = c(16, 16), col = c("blue", "red"), bty = "n", cex = 1.5)
dev.off()




if (T)
{
  
  
  tmpsp = seq(0, by = 0.05, len = 10)
  tmp = dgamma(tmpsp, shape = 2, scale = 0.2)
  tmp = tmp /sum(tmp)
  tmp = rev(cumsum(rev(tmp)))
  
  
  pdf("../figure/pointOutOfPmfScore.pdf", width = 10, height = 10 * 0.618)
  par(mar = c(4.2, 5, 0, 0), family = "serif")
  plot(x = tmpsp, y = tmp, type = "h", xlab = "Damage ratio", ylab = "P", 
       bty = "L", cex.lab = 2, cex.axis = 1.5, las = 1, ylim = c(0, max(tmp) * 1.45),
       xlim = c(0, 0.7))
  axis(side = 1, at = tmpsp[length(tmpsp)], 
       labels = latex2exp::TeX('$x_{\\max}$'), cex.axis = 1.5)
  text(x = tmpsp[length(tmpsp)], y = tmp[length(tmp)] * 2, labels = 
         latex2exp::TeX('$p_{\\max}$'), cex = 1.5)
  text(x = 0.075, y = 0.001, labels = latex2exp::TeX('$\\Delta x$'), cex = 1.5 )
  lines(x = 0.325, y = 0, type = "p", cex = 1.5, pch = 16, col = "blue")
  pmid = mean(tmp[7:8])
  lines(x = 0.6, y = 0, type = "p", cex = 1.5, pch = 16, col = "red")
  legend("topleft", legend = c("u: claim within PMF support", 
                               "v: claim out of PMF support",
                               "CCDF", 
                               "Exponential CCDF"
                               ),
         pch = c(16, 16, NA, NA), col = c("blue", "red", "black", "olivedrab3"), 
         bty = "n", cex = 1.5, lwd = c(NA, NA, 2, 2))
  # legend("right", legend = c("CCDF", "Exponential CCDF"),
  #        lwd = c(2, 2), col = c("black", "olivedrab3"), 
  #        bty = "n", cex = 1.5, lty = c('solid', 'dashed'))
  lda = -log(tmp[length(tmp)]) / tmpsp[length(tmpsp)]
  lines(x = seq(0, tmpsp[length(tmpsp)] * 2, len = 100),
        y = 1 - pexp(q = seq(0, tmpsp[length(tmpsp)] * 2, len = 100), rate = lda),
        col = "olivedrab3", lty = 'dashed', lwd = 2)
  lines(x = tmpsp[length(tmpsp)], y = tmp[length(tmp)], pch = 4, cex = 1.5,
        col = "darkgreen", type = "p", lwd = 2)
  dev.off()
  
  
}








































