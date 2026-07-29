suppressWarnings(library(ggplot2))

test_that("chromatogram_raw_plot returns expected data type", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  result <- chromatogram_raw_plot(
    raw.df = raw_builtin_file,
    dxf_path = builtin_file,
    title = "test",
    write_pdf = FALSE
  )

  expect_null(result)
})
