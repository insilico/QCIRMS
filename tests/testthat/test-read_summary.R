

test_that("function returns expected data type", {
  invisible(capture.output(result <- read_summary(builtin_file)))

  expect_s3_class(result, "summaryDefault")
})
