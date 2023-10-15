set.seed(42)
X = data.frame(val = sort(runif(30)), P = c(runif(14), 300, runif(15))) # Input 
# PMF's column names do not have to be "val" and "P".


X$P = X$P / sum(X$P)
plot(X, type = "h") # Made a close-to-degenerate PMF to increase the chance 
# of algorithm entering Stage II. See the publication referenced in Details.


newSupport = sort(c(X$val[1] - runif(1) * 0.05,
                    runif(20), X$val[nrow(X)] + runif(1) * 0.05))


fourPointRegridDetails = keyALGs::rglr4splitDetails(X, newSupport)


fourPointRegridDetails




