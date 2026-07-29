

test_that("function returns expected data type", {
  result <- vendor_info(builtin_file)

  expect_s3_class(result, "data.frame")
})
