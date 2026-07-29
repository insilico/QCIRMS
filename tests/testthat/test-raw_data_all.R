
test_that("function returns expected data type", {
  result <- raw_data_all(builtin_dxf_folder_vector_list)

  expect_type(result, "list")

})
