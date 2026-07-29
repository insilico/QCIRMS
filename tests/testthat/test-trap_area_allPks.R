
test_that("function returns expected data type", {

  result <- suppressWarnings(
    trap_area_allPks(raw.df = raw_builtin_file, vend.df = vend_builtin_file, mV.rawName = "v44.mV")
  )

  expect_s3_class(result, "data.frame")
})
