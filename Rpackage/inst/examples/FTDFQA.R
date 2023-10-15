# ==============================================================================
# `tmpParam` comes with the package. It is a matrix of 4 rows. Each column
# contains (a, b, c, lm1). `lm1` is the limited mean.
# ==============================================================================
data(tmpParam, package = "NGFMfitDistr")


# ==============================================================================
# `P0forTest` comes with the package.
# ==============================================================================
data(P0forTest, package = "NGFMfitDistr")
P0 = P0forTest


# ==============================================================================
# Solve the scale parameter d.
# ==============================================================================
d = solve_d(tmpParam, eps = 1e-8, maxit = 100)
abcd = rbind(tmpParam[1:3, , drop = FALSE], d)
dimnames(abcd) = NULL


# ==============================================================================
# Let tentative max be equal to the 99-th percentile --- if this percentile
# is below 1.
# ==============================================================================
tailProbThreshold = 1 - 0.99


# ==============================================================================
# We do not want the max to be too low.
# ==============================================================================
maxLowerBound = 0.05


mdrs = c(1:1e4L, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5


# ==============================================================================
# Compute the maxes.
# ==============================================================================
maxes = pmin(1, actuar_qtrbeta(1 - tailProbThreshold / (1 - P0), abcd))
maxes[!is.finite(maxes)] = 0 # In case of errors.
maxes = pmax(maxLowerBound, maxes)
ind = longestNonDecreasingSubseq(maxes)
maxFun = splinefun(x = mdrs[ind], y = maxes[ind], method = "hyman")
maxes = maxFun(mdrs)


# ==============================================================================
# Create TrB parameter table.
# ==============================================================================
TrBtable = rbind(mdrs, maxes, P0, abcd)


# ==============================================================================
# Discretization
# ==============================================================================
pmftablelist = FTDFQA(
  TrBtable, 
  supportSize = c(42L, 64L), 
  regridMethod = "lr",
  RIBlib = "Numerical Recipes",
  outputProbsInRows = TRUE,
  fineDiscretizationSize = 2000,
  maxCore = 1000, 
  verbose = TRUE, 
  downScaleSupport = FALSE)





