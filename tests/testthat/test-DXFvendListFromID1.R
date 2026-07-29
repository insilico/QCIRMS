
test_that("function returns expected data type", {
  builtin_vend_df <- do.call(
    rbind,
    combineVendFileInfo(
      path = builtin_dxf_folder,
      outputID = "test",
      outPath = tempdir()
    )
  )

  result <- QCIRMS:::DXFvendListFromID1(
    all.df = builtin_vend_df,
    standName.vec = unique(builtin_vend_df$Identifier1)
  )

  expect_type(result, "list")
})

