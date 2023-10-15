X = data.frame(val = sort(runif(100)), P = runif(100))
X$P = X$P / sum(X$P)
Y = data.frame(val = sort(runif(200)), P = runif(200))
Y$P = Y$P / sum(Y$P)
XplusY = keyALGs::como(X, Y, rgMethod = "") # Do not regrid.


corrImposed = (keyALGs::Var(XplusY) - (keyALGs::Var(X) + keyALGs::Var(Y))) / 
  (2 * sqrt(keyALGs::Var(X) * keyALGs::Var(Y)))


cat("Result PMF has only nonnegative probabilities = ",
    min(XplusY$P) >= 0, "\n",
    "Sum of result PMF's probabilities = ",
    sum(XplusY$P), "\n",
    "Theoretical mean of the result PMF = ",
    keyALGs::Mean(X) + keyALGs::Mean(Y), "\n",
    "Actual mean of the result PMF = ",
    keyALGs::Mean(XplusY), "\n",
    "Maximum correlation that can be imposed between X and Y = ",
    corrImposed, "\n", sep = "")


XplusYregrid = keyALGs::como(X, Y, N = 50, rgMethod = "r4")
plot(XplusYregrid, type = "h")
