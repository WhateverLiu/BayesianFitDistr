# X's support must be regular except for the last step which can be less.
X = data.frame(val = seq(runif(1), runif(1) + 1, len = 100), P = runif(100))
stepSize = X$val[2] - X$val[1]
X = rbind(X, c( X$val[nrow(X)] + stepSize * runif(1), runif(1) ) )
X$P = X$P / sum(X$P)


# Y's support must be regular except for the last step which can be less.
Y = data.frame(val = seq(runif(1), runif(1) + 1, len = 200), P = runif(200))
stepSize = Y$val[2] - Y$val[1]
Y = rbind(Y, c( Y$val[nrow(Y)] + stepSize * runif(1), runif(1) ) )
Y$P = Y$P / sum(Y$P)


XplusY = keyALGs::convSplitAtom4productsIrregular(
  X, Y, regridMethod = "r4", N = 256)


cat("Result PMF has only nonnegative probabilities = ",
    min(XplusY$P) >= 0, "\n",
    "Sum of result PMF's probabilities = ",
    sum(XplusY$P), "\n",
    "Theoretical mean of the result PMF = ",
    keyALGs::Mean(X) + keyALGs::Mean(Y), "\n",
    "Actual mean of the result PMF = ",
    keyALGs::Mean(XplusY), "\n",
    "Theoretical variance of the result PMF = ",
    keyALGs::Var(X) + keyALGs::Var(Y), "\n",
    "Actual variance of the result PMF = ",
    keyALGs::Var(XplusY), "\n", 
    "Unique step sizes on the result PMF's support = ",
    paste0(unique(diff(XplusY$val)), collapse = ", "), ". Here, ",
    "the unique step sizes could be more than two in the printout, ",
    "but it is due to numeric precision limit.\n",
    sep = "")


par(mfrow = c(1, 3))
ylim = c(0, max(max(X$P), max(Y$P), max(XplusY$P)))
plot(X, type = "h", bty = "L", ylim = ylim, xlab = "", ylab = "")
plot(Y, type = "h", bty = "L", ylim = ylim, xlab = "", ylab = "")
plot(XplusY, type = "h", bty = "L", xlab = "", ylab = "", ylim = ylim)
