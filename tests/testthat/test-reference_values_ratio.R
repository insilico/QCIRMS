

test_that("function returns expected data type", {
  result <- reference_values_ratio(c(builtin_file))

  expect_s3_class(result, "data.frame")
})

