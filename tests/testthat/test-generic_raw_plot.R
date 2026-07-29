test_that("function returns expected data type", {
  pdf(NULL) #prevents a plot from being stored 
  on.exit(dev.off(), add = TRUE)

  result <- generic_raw_plot(raw.df = raw_builtin_file, title = "test", write_pdf = FALSE)

  expect_null(result)
})
