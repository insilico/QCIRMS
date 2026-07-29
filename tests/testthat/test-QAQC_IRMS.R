

test_that("no data passes checks with builtin_dxf", {
  
  invisible(capture.output(result <- QCIRMS:::QAQC_IRMS(
    unfilteredPath = builtin_dxf_folder,
    expRef.df = test_exp_df,
    diff.t = 25,
    checkIntStand = FALSE,
    standAcceptedVals_vec = c(1, 2, 3),
    standAcceptedSD_vec = c(0.2, 0.2, 0.2),
    dataName = "test",
    maxPkNum = 25,
    expectedNonSampPks = 1,
    sdCrefIso.thresh = 100,
    sdOrefIso.thresh = 100,
    checkRelDiffIntensity = FALSE,
    refSamp_relDiffInt.thresh = 1,
    ref_relDiffInt.thresh = 1,
    sdCsampIso.thresh = 100,
    sdOsampIso.thresh = 100,
    filter_flushPk = FALSE,
    flushExpT = 135,
    flushTint = 15,
    filter_first_sampPk = FALSE,
    firstSampExpT = 275,
    firstSampTint = 15,
    firstSampExpPkNr = 2,
    outPath = tempdir(),
    verbose = FALSE
    )))

  
  expect_null(result)
})

#add test for .dxf files that pass qaqc check

