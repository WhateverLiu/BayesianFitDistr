library(NGFMfitDistr)
data(NewZealandCvgA, package = "NGFMfitDistr")
str(NewZealandCvgA)


mdrs = c(1:10000, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5
dat = data.frame(MDR = NewZealandCvgA$MDR, CDR = NewZealandCvgA$CDR)
dat = dat[sample(nrow(dat)), , drop = FALSE]
dat = dat[order(dat$MDR), , drop = FALSE]


p0ens = ensembleP0(
  # ============================================================================
  # Users need to ensure data rows have been ordered by MDR and CDR.
  mdr = NewZealandCvgA$MDR,
  cdr = NewZealandCvgA$CDR,
  # ============================================================================
  windowSize = round(nrow(NewZealandCvgA) * 0.01),
  slidingSpeed = round(nrow(NewZealandCvgA) * 1e-4),
  sampleSize = round(nrow(NewZealandCvgA) * (1 - 1 / exp(1))),
  NsampleSets = 100L,
  targetMDRs = mdrs,
  interpolationMethod = "hyman",
  randomSeed = 42L,
  maxCore = parallel::detectCores(),
  tempDir = "../tempFiles/CharlieTempMP/C"
)


# Plot the ensemble average against MDR.
plot(x = mdrs, y = rowMeans(as.data.frame(p0ens)), type = "l", 
     xlab = "MDR", ylab = expression(P[0]), bty = "L", las = 1)
