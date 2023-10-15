




pdf("../figure/noExtrapolationBiasInScore.pdf", width = 10, height = 10 * 0.618)
x = seq(0, 12, len = 100)
p = dgamma(x, shape = 9, scale = 0.5)
d1 = data.frame(val = x, P = p)
d2 = data.frame(val = seq(5, 6, len = 100), P = p)
d1$P = d1$P / diff(d1$val[1:2])
d2$P = d2$P / diff(d2$val[1:2])
xlim = c(-1, 13)
ylim = range(c(0, max(d1$P), max(d2$P)))
par(mar = c(4.2, 5, 0, 0), family = "serif")
plot(d1, type = 'l', xlim = xlim, ylim = ylim, bty = "L", xlab = "Value", 
     ylab = "Density", las = 1, cex.lab = 2, cex.axis = 1.5, lwd = 2)
lines(d2, col = "blue", lwd = 2)
lines(x = -0.5, y = 0, type = "p", col = "red", pch = 19, cex = 1)
lines(x = 12.5, y = 0, type = "p", col = "red", pch = 19, cex = 1)
legend("topright", legend = c(expression(f[1]), expression(f[2]), 
                              expression("Realization")), cex = 2,
       bty = "n", col = c("black", "blue", "red"), lwd = c(2, 2, NA),
       pch = c(NA, NA, 19))
dev.off()






















