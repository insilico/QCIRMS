

test_that("function returns expected data type", {
  result <- raw_data(builtin_file)

  expect_s3_class(result, "data.frame")
})
