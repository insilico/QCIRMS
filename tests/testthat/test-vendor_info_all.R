

test_that("function returns expected data type", {
  result <- vendor_info_all(c(builtin_file))
  expect_type(result, "list")
  
})

