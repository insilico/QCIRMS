
# combList<-combineVendFileInfo(builtin_dxf_folder)
# 
# processed_dxf<-QCIRMS:::removeFailedAnalysesDXF(sepList=combList,
#                                          expRef.df=test_exp_df,
#                                          maxPkNum=18,
#                                          expRefPkNr=dim(test_exp_df)[1],
#                                          diff.t=10,
#                                          sdCrefIso.thresh=0.1, 
#                                          expectedNonSampPks=7,
#                                          relDiffInt.thresh=0.1, 
#                                          amplName="Ampl44",
#                                          sdOrefIso.thresh=0.1,
#                                          sdCsampIso.thresh=0.3,
#                                          sdOsampIso.thresh=0.2,
#                                          filter_flushPk = T,
#                                          filter_first_sampPk = T,
#                                          flushExpT = 135,
#                                          flushTint = 15,
#                                          firstSampExpT = 275,
#                                          firstSampTint = 15,
#                                          firstSampExpPkNr = 42,
#                                          verbose=F)
# 


processed_dxf <- vector("list", 8)
processed_dxf[[5]] <- data.frame(failed_sampProc_Analysis = 99)
processed_dxf[[6]] <- data.frame(failed_sampIso_Analysis = 99)
processed_dxf[[7]] <- list(data.frame(Analysis = 1, PeakNr = 1))
processed_dxf[[8]] <- list(data.frame(Analysis = 1, PeakNr = 2))

refsList <- processed_dxf[[7]]
sampsList <- processed_dxf[[8]]

test_that("function returns expected data type", {
  #refsList <- list(data.frame(Analysis = 1, PeakNr = 1))
  #sampsList <- list(data.frame(Analysis = 1, PeakNr = 2))

  result <- QCIRMS:::analysisNumsEqualDXF(refsList = refsList, sampsList = sampsList)
  expect_type(result, "logical")
})
