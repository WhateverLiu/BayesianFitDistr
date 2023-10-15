

intensities_2010_09_04 = readxl::read_xlsx(
  "../data/preProcessData/Updated/ClaimsWithIntensity_2010-09-04_UpdatedwRSMDR.xlsx")
intensities_2011_02_22 = readxl::read_xlsx(
  "../data/preProcessData/Updated/ClaimsWithIntensity_2011-02-22_UpdatedwRSMDR.xlsx")


claimsData <- rbind(intensities_2010_09_04, intensities_2011_02_22)


names(claimsData)[match("RS.MDR.A", names(claimsData))] <- "GUMMDR_A"
names(claimsData)[match("MDR.A", names(claimsData))]  <- "GUCMDR_A"
names(claimsData)[match("RS.MDR.C", names(claimsData))] <- "GUMMDR_C"
names(claimsData)[match("MDR.C", names(claimsData))]  <- "GUCMDR_C"
names(claimsData)[match("Lat...28", names(claimsData))] <- "GeoLat"
names(claimsData)[match("Lon...29", names(claimsData))] <- "GeoLong"
claimsData = as.data.frame(claimsData)
claimsData = claimsData[c(
  "ClaimID","EventDate", "GeoLat","GeoLong",
  "OCC", "GUMMDR_A", "GUCMDR_A","GUMMDR_C","GUCMDR_C")]


claimsData$ClaimID = as.integer(claimsData$ClaimID)
claimsData$OCC = as.integer(claimsData$OCC)
claimsData$GUCMDR_C = as.numeric(claimsData$GUCMDR_C)
claimsData$GUCMDR_A = as.numeric(claimsData$GUCMDR_A)


data.table::fwrite(claimsData, file = "../data/cleanedClaimsData.csv")


claimsData = claimsData[claimsData$OCC %in% as.integer(c(
  "301", "302","303", "304", "305", "306", "307")), ]


cvgAdata = claimsData[, 1:7]
cvgAdata = cvgAdata[
  is.finite(cvgAdata$GUMMDR_A) & is.finite(cvgAdata$GUCMDR_A), ]
cvgAdata = cvgAdata[cvgAdata$GUMMDR_A <= 1 & cvgAdata$GUMMDR_A > 0 &
  cvgAdata$GUCMDR_A <= 1 & cvgAdata$GUCMDR_A >= 0, ]


cvgCdata = claimsData[, c(1:5, 8, 9)]
cvgCdata = cvgCdata[
  is.finite(cvgCdata$GUMMDR_C) & is.finite(cvgCdata$GUCMDR_C), ]
cvgCdata = cvgCdata[cvgCdata$GUMMDR_C <= 1 & cvgCdata$GUMMDR_C > 0 &
                      cvgCdata$GUCMDR_C <= 1 & cvgCdata$GUCMDR_C >= 0, ]


dir.create('../data', showWarnings = F)


cvgData = cvgAdata
cvgData$MDR = cvgData$GUMMDR_A
cvgData$CDR = cvgData$GUCMDR_A
cvgData$GUMMDR_A = NULL
cvgData$GUCMDR_A = NULL
save(cvgData, file = "../data/cvgAdata.Rdata")


cvgData = cvgCdata
cvgData$MDR = cvgData$GUMMDR_C
cvgData$CDR = cvgData$GUCMDR_C
cvgData$GUMMDR_C = NULL
cvgData$GUCMDR_C = NULL
save(cvgData, file = "../data/cvgCdata.Rdata")



















