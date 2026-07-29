
test_that("function returns expected data type", {

  oldwd <- getwd()
  setwd(tempdir())
  on.exit(setwd(oldwd), add = TRUE)

  invisible(capture.output(builtin_vend_df <- do.call(rbind, combineVendFileInfo(
      path = builtin_dxf_folder,
      outputID = "test",
      outPath = tempdir()
    )
    )))


  invisible(capture.output(result <- internal_standards_summary(
    data.df = builtin_vend_df,
    dataName = "test",
    standName.vec = c("0.1 m Na2SO4", "2.0 m Na2SO4", "5.0 m NaCl"),
    standAcceptedVals_vec = c(-5.75, -6.10, -5.80),
    standAcceptedSD_vec = c(2, 2, 2)
  )))
  
  expect_type(result, "list")
  expect_length(result, 4)
  expect_equal(result[[2]]$standard, c("0.1 m Na2SO4", "2.0 m Na2SO4", "5.0 m NaCl"))
})
