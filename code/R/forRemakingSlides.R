


exampleDistr = as.data.frame(empDistrsLists[[1]][[500]])
exampleDistrUpscale = as.data.frame(empDistrsLists[[1]][[15000]])
save(exampleDistr, exampleDistrUpscale, file = "../data/exampleDistr.Rdata")


# ==============================================================================
# Show how conditional PMF is made.
# ==============================================================================
if (T)
{
  load("../data/exampleDistr.Rdata")
  pdf("../figure/conditionalPMFexample.pdf", width = 10, height = 10 * 0.618)
  par(mar = c(4.2, 5, 0, 0), family = "serif", mfrow = c(2, 1))
  ylim = c(0, max(exampleDistr$P))
  plot(exampleDistr[-1,], type = "h", col = "blue", 
       ylim = c(0, max(exampleDistr$P)), bty = "L", cex.lab = 2, 
       cex.axis = 1.5, xlab = "", ylab = "P", las = 1)
  lines(x = c(0, 0), y = c(0, exampleDistr$P[1]), col = "red")
  coor = keyALGs::plotCoor(0.5, 0.9)
  text(x = coor[1], y = coor[2], labels = "Full empirical", cex = 1.8)
  coor = keyALGs::plotCoor(0.5, 0.7)
  text(x = coor[1], y = coor[2], labels = "MDR = 0.005", cex = 1.8)
  coor = keyALGs::plotCoor(0.5, 0.5)
  text(x = coor[1], y = coor[2], labels = "Mean = 0.054", cex = 1.8)
  coor = keyALGs::plotCoor(0.5, 0.3)
  text(x = coor[1], y = coor[2], labels = expression(P[0](MDR==0.005)==0.3033), cex = 1.8)
  tmp = exampleDistr[-1,]
  tmp$P = tmp$P / sum(tmp$P)
  plot(tmp, ylim = ylim, xlim = c(0, max(tmp$val)), type = "h", bty = "L",
       col = "blue", cex.lab = 2, cex.axis = 1.5, xlab = "Damage ratio", 
       ylab = "", las = 1)
  coor = keyALGs::plotCoor(0.5, 0.9)
  text(x = coor[1], y = coor[2], 
       labels = "Conditional empirical | Damage ratio > 0", cex = 1.8)
  coor = keyALGs::plotCoor(0.5, 0.7)
  text(x = coor[1], y = coor[2], labels = "MDR = 0.005", 
       cex = 1.8)
  coor = keyALGs::plotCoor(0.5, 0.5)
  text(x = coor[1], y = coor[2], labels = "Target mean for TrB = 0.005 / (1 - 0.3033)", 
       cex = 1.8)
  coor = keyALGs::plotCoor(0.5, 0.3)
  text(x = coor[1], y = coor[2], labels = "Mean = 0.054 / (1 - 0.3033)", 
       cex = 1.8)
  dev.off()
  
}








# ==============================================================================
# How the scaling of PMFs, downscale
# ==============================================================================
if (T)
{
  
  load("../data/exampleDistr.Rdata")
  tmp = exampleDistr[-1,]
  mx = round(max(tmp$val), 4)
  tmp$P = tmp$P / sum(tmp$P)
  exampleDistrMDR = 5e-3
  ylim = c(0, max(exampleDistr$P))
  pdf("../figure/downScaleExample.pdf", width = 10, height = 10 * 0.618)
  par(mar = c(4.2, 5, 0, 0), family = "serif")
  for (mdr in c(5e-3, 4e-3, 3e-3, 2e-3, 1e-3))
  {
    tmp = exampleDistr[-1,]
    tmp$P = tmp$P / sum(tmp$P)
    tmp$val = tmp$val * (mdr / exampleDistrMDR)
    plot(tmp, ylim = ylim, type = "h", bty = "L",
         col = "blue", cex.lab = 2, cex.axis = 1.5, xlab = "Damage ratio", 
         ylab = "P", las = 1)
    coor = keyALGs::plotCoor(0.5, 0.9)
    text(x = coor[1], y = coor[2], 
         labels = "Conditional empirical | Damage ratio > 0", cex = 1.8)
    coor = keyALGs::plotCoor(0.5, 0.8)
    text(x = coor[1], y = coor[2], labels = paste0("MDR = ", mdr), 
         cex = 1.8)
    coor = keyALGs::plotCoor(0.5, 0.7)
    if ( mdr == 5e-3 ) lb = paste0("Max = ", mx)
    else lb = paste0("Max = ", mx, " \u00D7 (", mdr, " / 0.005)")
    
    text(x = coor[1], y = coor[2], labels = lb, cex = 1.8)
    expr = eval(parse(text = paste0(
      "expression('Target mean for TrB = ", mdr, " / [1 -'~P[0]~'('~",
      mdr, "~')]')")))
    coor = keyALGs::plotCoor(0.5, 0.6)
    text(x = coor[1], y = coor[2], labels = expr, cex = 1.8)
  }
  dev.off()
  
  
}




# ==============================================================================
# How the scaling of PMFs, upscale
# ==============================================================================
if (T)
{
  
  
  load("../data/exampleDistr.Rdata")
  exampleDistr = exampleDistrUpscale
  
  
  tmp = exampleDistr[-1,]
  mx = round(max(tmp$val), 4)
  tmp$P = tmp$P / sum(tmp$P)
  exampleDistrMDR = 0.5
  ylim = c(0, max(exampleDistr$P))
  pdf("../figure/upScaleExample.pdf", width = 10, height = 10 * 0.618)
  par(mar = c(4.2, 5, 0, 0), family = "serif")
  dlist = list()
  for (mdr in c(0.5, 0.6, 0.7, 0.8, 0.9))
  {
    tmp = exampleDistr[-1,]
    tmp$P = tmp$P / sum(tmp$P)
    tmp$val = tmp$val * (mdr / exampleDistrMDR)
    pmax = sum(tmp$P[tmp$val >= 1])
    tmp = tmp[tmp$val < 1, ]
    dlist[[length(dlist) + 1]] = data.frame(
      val = c(tmp$val, 1), P = c(tmp$P, pmax) )
  }
  ylim = c(0, max(unlist(lapply(dlist, function(x) max(x$P)))))
  i = 1L
  for (mdr in c(0.5, 0.6, 0.7, 0.8, 0.9))
  {
    tmp = dlist[[i]]
    plot(tmp, ylim = ylim, type = "h", bty = "L",
         col = "blue", cex.lab = 2, cex.axis = 1.5, xlab = "Damage ratio", 
         ylab = "P", las = 1, xlim = c(0, 1))
    coor = keyALGs::plotCoor(0.5, 0.9)
    text(x = coor[1], y = coor[2], 
         labels = "Conditional empirical | Damage ratio > 0", cex = 1.8)
    coor = keyALGs::plotCoor(0.5, 0.8)
    text(x = coor[1], y = coor[2], labels = paste0("MDR = ", mdr), 
         cex = 1.8)
    coor = keyALGs::plotCoor(0.5, 0.7)
    if ( mdr == 0.5 ) lb = paste0("Max = ", mx)
    else lb = paste0("Max = min[1, ", mx, " \u00D7 (", mdr, " / 0.5)]")
    text(x = coor[1], y = coor[2], labels = lb, cex = 1.8)
    expr = eval(parse(text = paste0(
      "expression('Target mean for TrB = ", mdr, " / [1 -'~P[0]~'('~",
      mdr, "~')]')")))
    coor = keyALGs::plotCoor(0.5, 0.6)
    text(x = coor[1], y = coor[2], labels = expr, cex = 1.8)
    i = i + 1L
  }
  dev.off()
  
  
}




source("../code/R/rfuns.R")
# ==============================================================================
# Adapt franchise deductible
# ==============================================================================
if (T)
{
  
  Npoint = 100L- 1L
  ded = 0.1
  p0star = 0.5
  shape = 3; scale = 2
  set.seed(789)
  p = dgamma(seq(0, 20, len = Npoint), shape = shape, scale = scale)
  distr = data.frame(val = seq(0, 1, len = Npoint), P = p)
  truncatedDistr = distr[distr$val >= ded, ]
  truncatedDistrJagged = truncatedDistr
  truncatedDistrJagged$P = truncatedDistrJagged$P + 
    runif(nrow(truncatedDistrJagged), -1, 1) * 
    (max(truncatedDistrJagged$P) * 0.075)
  truncatedDistrJagged$P = pmax(0, truncatedDistrJagged$P)
  
  
  
  
  pdf("../figure/explainFranchiseDed.pdf", width = 10, height = 10 * 0.618)
  
  
  par(mar = c(4.2, 5, 0, 0), family = "serif", mfrow = c(2, 2))
  
  
  tmp = truncatedDistrJagged
  # tmp = data.frame(val = c(0, tmp$val), P = c(p0star, tmp$P))
  tmpHist = wrapAsHist(tmp)
  tmpHist$breaks = c(0, tmpHist$breaks)
  tmpHist$density = c(p0star, tmpHist$density)
  colr = scales::alpha("blue", 0.5)
  colr = c(scales::alpha("red", 0.5), rep(colr, length(tmpHist$density)))
  
  
  plot(tmpHist, freq = F, border = NA, col = colr, ylab = "P", xlab = "", 
       main = "", las = 1, cex.lab = 2, cex.axis = 1.5)
  tmp = paste0('$P^{*}_{0}= ', p0star, '$')
  legend("topright", legend = c(
    latex2exp::TeX(tmp), latex2exp::TeX('Main part')), 
         bty = "n" , cex = 2, col = colr[1:2], pch = 15)
  legend("center", legend = paste0("Franchise deductible = ", ded ), cex = 1.5, 
         bty = "n")
  legend("right", legend = "(1)", bty = "n", cex = 2)
  
  
  plot(tmpHist, freq = F, border = NA, col = colr, ylab = "", xlab = "", 
       main=  "", las = 1, cex.lab = 2, cex.axis = 1.5)
  lines(distr, col = "black", lwd = 2)
  legend("topright", legend = "Tentative TrB fit", bty = "n", cex = 2, lwd = 2, 
         col = "black")
  legend("right", legend = "(2)", bty = "n", cex = 2)
  
  
  plot(tmpHist, freq = F, border = NA, col = colr, ylab = "", 
       xlab = "Damage ratio", main=  "", las = 1, cex.lab = 2, cex.axis = 1.5)
  lines(distr, col = "black", lwd = 2)
  tmp = distr[1:max(which(distr$val < ded) + 1L), ]
  # tmp2 = tmpHist$breaks[2]
  # tmp = rbind( tmp, c(tmp2, tmp$P[nrow(tmp)]), c(tmp2, 0)   )
  tmp2 = max(tmp$val)
  tmp = rbind(tmp, c(tmp2, 0) )
  polygon(tmp, col = "black", border = NA)
  area = round(sum(tmp$P) * ded, 3)
  legend("topright", legend = paste0(
    "Area = ", area  ), bty = "n", pch = 15, col = "black", cex = 2)
  legend("right", legend = "(3)", bty = "n", cex = 2)
  
  
  colr[1] = "white"
  p0 = p0star - area
  tmpHist$density = tmpHist$density * (  (1 - p0) / (1 - p0star) )
  plot(tmpHist, freq = F, border = NA, col = colr, 
       ylab = "", xlab = "", main=  "", las = 1, cex.lab = 2, cex.axis = 1.5)
  lines(x = 0, y = p0star - area, type = "h", col = "red", lwd = 2 )
  legend("top", legend = "Final empirical to be fitted by TrB", bty= "n", cex = 2)
  tmp = paste0(
    '$P_0 \\leftarrow P^*_0 - ', area, '=', p0star - area, '$')
  legend("center", legend = latex2exp::TeX(tmp), bty= "n", cex = 1.5, 
         col = "red", lwd = 2)
  legend("right", legend = "(4)", bty = "n", cex = 2)
  
  
  
  dev.off()
  
  
  
  
  
  
  
}




# ==============================================================================
# Adapt ordinary deductible
# ==============================================================================
if (T)
{
  
  Npoint = 100L- 1L
  ded = 0.1
  p0star = 0.5
  shape = 3; scale = 2
  set.seed(789)
  p = dgamma(seq(0, 20, len = Npoint), shape = shape, scale = scale)
  distr = data.frame(val = seq(0, 1, len = Npoint), P = p)
  truncatedDistr = distr[distr$val >= ded, ]
  truncatedDistrJagged = truncatedDistr
  truncatedDistrJagged$P = truncatedDistrJagged$P + 
    runif(nrow(truncatedDistrJagged), -1, 1) * 
    (max(truncatedDistrJagged$P) * 0.075)
  truncatedDistrJagged$P = pmax(0, truncatedDistrJagged$P)
  
  
  pdf("../figure/explainOrdinaryDed.pdf", width = 10, height = 10 * 0.618)
  
  
  par(mar = c(4.2, 5, 0, 0), family = "serif")
  layout(rbind(c(1, 1, 2, 2),
               c(3, 3, 4, 4),
               c(5, 6, 6, 7)))
  
  
  tmp = truncatedDistrJagged
  tmpHist = wrapAsHist(tmp)
  tmpHist$breaks = c(0, tmpHist$breaks)
  tmpHist$density = c(p0star, tmpHist$density)
  colr = scales::alpha("blue", 0.5)
  colr = c(scales::alpha("red", 0.5), rep(colr, length(tmpHist$density)))
  
  
  tmp = tmpHist
  tmp$density = tmp$density[-1]
  tmp$breaks = tmp$breaks[-1] - ded
  plot(tmp, freq = F, border = NA, col = colr[-1], ylab = "P", xlab = "", 
       main = "", las = 1, cex.lab = 2, cex.axis = 1.5, ylim = c(0, p0star), 
       xlim = c(0, 1))
  lines(x = 0, y = p0star, type ="h", lwd = 2, col = "red")
  legend("topleft", legend = paste0("Ordinary\ndeductible = ", ded), bty = "n", 
         cex = 2)
  tmp = paste0('$P^{*}_{0}= ', p0star, '$')
  legend("topright", legend = c(
    latex2exp::TeX(tmp), latex2exp::TeX('Main part')), 
    bty = "n" , cex = 2, col = colr[1:2], pch = 15)
  legend("right", legend = "(1)", bty = "n", cex = 2)
  
  
  tmp = tmpHist
  tmp$density = tmp$density[-1]
  tmp$breaks = tmp$breaks[-1]
  plot(tmp, freq = F, border = NA, col = colr[-1], ylab = "", xlab = "", 
       main = "", las = 1, cex.lab = 2, cex.axis = 1.5, ylim = c(0, p0star), 
       xlim = c(0, 1))
  lines(x = 0, y = p0star, type ="h", lwd = 2, col = "red")
  legend("top", legend = paste0("Right shift by ", ded), bty = "n", cex = 2)
  legend("right", legend = "(2)", bty = "n", cex = 2)
  
  
  # plot(tmpHist, freq = F, border = NA, col = colr, ylab = "", xlab = "", 
  #      main = "", las = 1, cex.lab = 2, cex.axis = 1.5)
  # # tmp = paste0('$P^{*}_{0}= ', p0star, '$')
  # # legend("topright", legend = c(
  # #   latex2exp::TeX(tmp), latex2exp::TeX('Main part')), 
  # #   bty = "n" , cex = 2, col = colr[1:2], pch = 15)
  # legend("right", legend = "(3)", bty = "n", cex = 2)
  
  
  tmpcolr = colr; tmpcolr[1] = "white"
  plot(tmpHist, freq = F, border = NA, col = tmpcolr, ylab = "", xlab = "", 
       main = "", las = 1, cex.lab = 2, cex.axis = 1.5)
  lines(distr, col = "black", lwd = 2)
  lines(x = 0, y = p0star, type ="h", lwd = 2, col = "red")
  legend("topright", legend = "Tentative TrB fit", bty = "n", 
         cex = 2, lwd = 2, col = "black")
  legend("right", legend = "(3)", bty = "n", cex = 2)
  
  
  tmpcolr = colr; tmpcolr[1] = "white"
  plot(tmpHist, freq = F, border = NA, col = tmpcolr, ylab = "", 
       xlab = "", main=  "", las = 1, cex.lab = 2, cex.axis = 1.5)
  lines(distr, col = "black", lwd = 2)
  lines(x = 0, y = p0star, type ="h", lwd = 2, col = "red")
  tmp = distr[1:max(which(distr$val < ded) + 1L), ]
  tmp2 = max(tmp$val)
  tmp = rbind(tmp, c(tmp2, 0) )
  polygon(tmp, col = "black", border = NA)
  area = round(sum(tmp$P) * ded, 3)
  legend("topright", legend = paste0(
    "Area = ", area  ), bty = "n", pch = 15, col = "black", cex = 2)
  legend("right", legend = "(4)", bty = "n", cex = 2)
  
  
  plot.new()
  
  
  colr[1] = "white"
  p0 = p0star - area
  tmpHist$density = tmpHist$density * (  (1 - p0) / (1 - p0star) )
  plot(tmpHist, freq = F, border = NA, col = colr, 
       ylab = "", xlab = "Damage ratio", main=  "", 
       las = 1, cex.lab = 2, cex.axis = 1.5)
  lines(x = 0, y = p0star - area, type = "h", col = "red", lwd = 2 )
  legend("top", legend = "Final empirical to be fitted by TrB", bty= "n", cex = 2)
  tmp = paste0('$P_0 \\leftarrow P^*_0 - ', 
               area, '=', p0star - area, '$')
  legend("center", legend = latex2exp::TeX(tmp), bty= "n", cex = 1.75, 
         col = "red", lwd = 2)
  legend("right", legend = "(5)", bty = "n", cex = 2)
  
  
  
  dev.off()
  
  
  
  
  
  
  
}


































