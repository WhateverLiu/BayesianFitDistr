X = data.frame(val = sort(runif(100)), P = runif(100))
X$P = X$P / sum(X$P)
Y = data.frame(val = sort(runif(200)), P = runif(200))
Y$P = Y$P / sum(Y$P)
XplusY = keyALGs::convBt(X, Y, eps = 1e-10)


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
    keyALGs::Var(XplusY), "\n", sep = "")
