

test_that("function returns expected data type", {
  result <- QCIRMS:::maxNumPeaksDXF(
    vend.df = vend_builtin_file,
    maxExpectedPks = 10,
    expectedNonSampPks = 5,
    verbose = FALSE
  )

  expect_type(result, "list")
})
