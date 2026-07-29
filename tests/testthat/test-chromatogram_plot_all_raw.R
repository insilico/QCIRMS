suppressWarnings(library(ggplot2))

test_that("chromatogram_plot_all_raw returns expected data type", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  result <- chromatogram_plot_all_raw(
    raw.list = raw_builtin_dxf_folder,
    result_path = tempdir(),
    data_path = builtin_dxf_folder,
    write_pdf = FALSE
  )

  expect_type(result, "character")
  expect_length(result, 1)
})
