X = data.frame(val = sort(runif(100)), P = runif(100)) # Input PMF's column
# names do not have to be "val" and "P".
X$P = X$P / sum(X$P)
newSupport = sort(c(runif(28), X$val[1], X$val[nrow(X)]))
Y = keyALGs::rglr(X, newSupport)


# Plot the two PMFs.
ylim = c(0, max(max(X$P), max(Y$P)))
par(mfrow = c(1, 2))
plot(X, type = "h", bty = "L", ylim = ylim)
legend("top", legend = expression(p[X]), bty = "n")
plot(Y, type = "h", bty = "L", ylim = ylim, ylab = "", xlab = "")
legend("top", legend = expression(p[Y]), bty = "n")
