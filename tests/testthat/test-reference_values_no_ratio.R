
test_that("function returns expected data type", {
  result <- reference_values_no_ratio(builtin_dxf_folder_vector_list)

  expect_s3_class(result, "data.frame")
})

