






test_that("function returns expected data type", {
  result <- combineVendFileInfo(builtin_dxf_folder)


  expect_type(result, "list")
})
