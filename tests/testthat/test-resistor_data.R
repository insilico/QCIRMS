
test_that("function returns expected data type", {
  result <- resistor_data(c(builtin_file))

  expect_s3_class(result, "data.frame")
  
})
