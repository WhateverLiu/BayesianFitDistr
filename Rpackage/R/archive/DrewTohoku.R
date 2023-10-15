TohokuReduced = data.table::setDF(data.table::fread("../tempFiles/TohokuReduced.csv"))
load('../tempFiles/allData2.RData')


TohokuReduced$V1 = NULL
TohokuReduced$X = NULL
TohokuReduced = TohokuReduced[
  TohokuReduced$Model_DRb_shk + TohokuReduced$Model_DRb_liq > TohokuReduced$Model_DRb_tsu, , drop = F]
ind = TohokuReduced$peril_of_damage %in% c(81L, 42L, 82L)
TohokuReduced = TohokuReduced[ind, , drop = F]
TohokuReduced = TohokuReduced[
  TohokuReduced$Model_DRb_shk + TohokuReduced$Model_DRc_shk  + 
    TohokuReduced$Model_DRb_liq > 0, , drop = F]
TohokuReduced = TohokuReduced[TohokuReduced$Distance_to_Fukushima > 30, , drop = F]
TohokuReduced = TohokuReduced[!(
  TohokuReduced$Limit_A_JPY < 1e-10 & TohokuReduced$Limit_C_JPY < 1e-10), , drop = F]
TohokuReduced$GUCDR_A = TohokuReduced$Payment_A_JPY / (
  TohokuReduced$SI_A_JPY * TohokuReduced$RI)
TohokuReduced$GUCDR_C = TohokuReduced$Payment_C_JPY / (
  TohokuReduced$SI_C_JPY * TohokuReduced$RI)
TohokuReduced = TohokuReduced[
  !(is.nan(TohokuReduced$GUCDR_A) & is.nan(TohokuReduced$GUCDR_C)), , drop = F]


cvgA = data.frame(MDR = TohokuReduced$Model_DRb_shk + TohokuReduced$Model_DRb_liq, 
                  CDR = TohokuReduced$GUCDR_A)
cvgA = cvgA[is.finite(cvgA$MDR) & is.finite(cvgA$CDR), , drop = F]
cvgA = cvgA[cvgA$MDR > 0 & cvgA$MDR < 1 & cvgA$CDR >= 0 & cvgA$CDR <= 1, , drop = F]
cvgA$CDR[cvgA$CDR < 0.05] = 0


cvgC = data.frame(MDR = TohokuReduced$Model_DRc_shk, CDR = TohokuReduced$GUCDR_C)
cvgC = cvgC[is.finite(cvgC$MDR) & is.finite(cvgC$CDR), , drop = F]
cvgC = cvgC[cvgC$MDR > 0 & cvgC$MDR < 1 & cvgC$CDR >= 0 & cvgC$CDR <= 1, , drop = F]
cvgC$CDR[cvgC$CDR < 0.05] = 0


save(cvgA, file = "../data/tohokuCvgA.Rdata")
save(cvgC, file = "../data/tohokuCvgC.Rdata")

