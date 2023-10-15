X = data.frame(val = seq(0, 16, len = 1000), P = dgamma(
  seq(0, 16, len = 1000), shape = 2, scale = 2)) # Input PMF's column 
# names do not have to be "val" and "P".


X$P = X$P / sum(X$P)
newGrid = seq(0, 16, len = 64)
Xnew = keyALGs::LMMdiscreteToDiscrete(X, newGrid)


# Means:
keyALGs::Mean(X); keyALGs::Mean(Xnew)
# 3.957082
# 3.957082


# Variances:
keyALGs::Var(X); keyALGs::Var(Xnew)
# 7.396614
# 7.396614
