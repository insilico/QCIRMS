checkRelDiffIntensity <- T

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

# adjust reference data for any samples removed
# list of two lists: one for ref peak dfs, one for samps
currFiltered<-QCIRMS:::rmRefDatDXF(currProc)
# want same number of analyses in both
lengthEqual<-length(currFiltered[[1]])==length(currFiltered[[2]])
#lengthEqual

if(length(currFiltered)>0){
  # check relative differences bt samp and ref intensity
  if(checkRelDiffIntensity && lengthEqual){

  test_that("function returns expected data type", {
  # ref.df <- data.frame(Analysis = 1, Identifier1 = "sample", PeakNr = c(1, 3), Ampl44 = c(100, 100))
  # samp.df <- data.frame(Analysis = 1, Identifier1 = "sample", PeakNr = c(2, 4), Ampl44 = c(101, 101))
  # currFiltered <- list(ref.df, samp.df)
    invisible(capture.output(result<-ref_samp_intensity_check(currFiltered=currFiltered,
                                         lengthEqual=lengthEqual,
                                         dataName="test",
                                         amplName="Ampl44",
                                         relDiffInt.thresh=0.2,
                                         outPath = tempdir(),
                                         verbose=F)))

  
  # result <- ref_samp_intensity_check(
  #   currFiltered = vend_builtin_file,
  #   lengthEqual = TRUE,
  #   dataName = "test",
  #   amplName = "Ampl44",
  #   relDiffInt.thresh = 0.2,
  #   outPath = tempdir(),
  #   verbose = FALSE
  # )

  expect_type(result, "list")
})
    
  }
}
