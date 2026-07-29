
test_that("function returns expected data type", {
  ref.list <- list(data.frame(Analysis = 1, PeakNr = c(1, 2), Ampl44 = c(100, 101)))
  sample.list <- list(data.frame(Analysis = 1, PeakNr = c(3, 4), Ampl44 = c(90, 91)))

  result <- QCIRMS:::threadRefSampListsToDF(ref.list = ref.list, sample.list = sample.list)

  expect_type(result, "list")
})

