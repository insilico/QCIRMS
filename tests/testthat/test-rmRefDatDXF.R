

combList<-combineVendFileInfo(builtin_dxf_folder, outputID = "test")

currProc <- QCIRMS:::removeFailedAnalysesDXF(sepList=combList,
                                                 expRef.df=test_exp_df,
                                                 maxPkNum=18,
                                                 expRefPkNr=dim(test_exp_df)[1],
                                                 diff.t=10,
                                                 sdCrefIso.thresh=0.1, 
                                                 expectedNonSampPks=7,
                                                 relDiffInt.thresh=0.1, 
                                                 amplName="Ampl44",
                                                 sdOrefIso.thresh=0.1,
                                                 sdCsampIso.thresh=0.3,
                                                 sdOsampIso.thresh=0.2,
                                                 filter_flushPk = F, #must be false
                                                 filter_first_sampPk = F, #must be false
                                                 flushExpT = 135,
                                                 flushTint = 15,
                                                 firstSampExpT = 275,
                                                 firstSampTint = 15,
                                                 firstSampExpPkNr = 42,
                                                 verbose=F)


test_that("function returns expected data type", {
  result <- QCIRMS:::rmRefDatDXF(currProc)

  expect_type(result, "list")
})
