
test_that("function returns expected data type", {
  result <- file_info(c(builtin_file))

  expect_s3_class(result, "data.frame")
})
