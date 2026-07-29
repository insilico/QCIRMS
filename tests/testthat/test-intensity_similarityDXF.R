
test_that("function returns expected data type", {
  result <- intensity_similarityDXF(
    vendAmpl = c(100, 101, 102),
    amplName = "Ampl44",
    peakNr.vec = 1:3,
    relDiffInt.thresh = 0.1,
    verbose = FALSE
  )

  expect_type(result, "list")
})
