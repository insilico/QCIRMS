

test_that("function returns expected data type", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  result <- generic_plot_all_raw(raw_list = raw_builtin_dxf_folder, write_pdf = FALSE)

  expect_null(result)
})
