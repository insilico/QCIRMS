
#change folder path

test_that("function returns expected data type", {

  result <- all_dxf_files(builtin_dxf_folder)

  expect_type(result, "character")
})

