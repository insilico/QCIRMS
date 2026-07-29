
test_that("function returns expected data type", {
  result <- isoR_similarityDXF(
    vend.df = vend_builtin_file,
    peakNr.vec = 1:5,
    sdC.thresh = 0.1,
    sdO.thresh = 0.1,
    verbose = FALSE
  )

  expect_type(result, "list")
})

