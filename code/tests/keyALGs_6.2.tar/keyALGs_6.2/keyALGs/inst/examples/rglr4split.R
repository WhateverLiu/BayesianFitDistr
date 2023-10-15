X = data.frame(val = sort(runif(100)), P = runif(100))  # Input PMF's column 
# names do not have to be "val" and "P".
X$P = X$P / sum(X$P)
newSupport = sort(c(runif(28), X$val[1], X$val[nrow(X)]))
linearRegridedX = keyALGs::rglr(X, newSupport)
fourPointRegridedX = keyALGs::rglr4split(X, newSupport)


# Linear regrid inflated variance. Four-point regrid preserved variance. 
cat("Mean(X), Var(X) = ", 
    keyALGs::Mean(X), ", ", keyALGs::Var(X), "\n",
    "Mean(linearRegridedX), Var(linearRegridedX) = ", 
    keyALGs::Mean(linearRegridedX), ", ", keyALGs::Var(linearRegridedX), "\n",
    "Mean(fourPointRegridedX), Var(fourPointRegridedX) = ", 
    keyALGs::Mean(fourPointRegridedX), ", ", keyALGs::Var(fourPointRegridedX), 
    sep = "")
