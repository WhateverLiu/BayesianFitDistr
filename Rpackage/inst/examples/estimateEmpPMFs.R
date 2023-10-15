library(NGFMfitDistr)
data(NewZealandCvgA, package = "NGFMfitDistr"); str(NewZealandCvgA)


# ==============================================================================
# A P0 model just for example.
# ==============================================================================
data(P0forTest, package = "NGFMfitDistr")


# ==============================================================================
# Make sliding windows.
# ==============================================================================
windows = computeWindows(
  xlen = nrow(NewZealandCvgA), 
  windowSize = round(nrow(NewZealandCvgA) * 0.01),
  start = 1L,
  speed = round(nrow(NewZealandCvgA) * 1e-4),
  returnList = TRUE
  )


# ==============================================================================
# Because NewZealandCvgA has been ordered by MDR and CDR, we need to reshuffle
# the rows of NewZealandCvgA and resort it by MDR only. This is because there
# are many different CDRs associated to the same MDR in this dataset.
# ==============================================================================
set.seed(123)
NewZealandCvgA = NewZealandCvgA[sample(nrow(NewZealandCvgA)), , drop = FALSE]
NewZealandCvgA = NewZealandCvgA[order(NewZealandCvgA$MDR), , drop = FALSE]
  

rst = estimateEmpPMFs(
    mdr = NewZealandCvgA$MDR, 
    cdr = NewZealandCvgA$CDR, 
    windows = windows, 
    P0 = P0forTest,
    targetMDRs = c(1:1e4L, seq(1e4L + 10L, 1e5L - 10L, by = 10L)) / 1e5,
    regridMethod = "lmm",
    nonDegenerateOldDistrs = NULL,
    sampleWeightOnOldDistrs = 0.5,
    empDistrSupportSize = 64L,
    maxCore = parallel::detectCores())


