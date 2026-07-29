
test_that("function returns expected data type; detects correct num of ref peaks
          detects wrong num of ref peaks", {

  invisible(capture.output(result <- QCIRMS:::reference_times_checkDXF(
    vend.df = vend_builtin_file,
    expectedPeak.num = nrow(test_exp_df),
    diff.t = 10,
    expectedStart = test_exp_df$Expected_Start,
    expectedRt = test_exp_df$Expected_Rt,
    expectedEnd = test_exp_df$Expected_End,
    verbose = TRUE
  )))

  expect_type(result, "list")
  expect_true(result[[1]])
  expect_equal(result[[2]]$PkNr, c(1, 2, 4, 5, 17))
  
  invisible(capture.output(fail_result <- QCIRMS:::reference_times_checkDXF(
    vend.df = vend_builtin_file,
    expectedPeak.num = 999,
    diff.t = 10,
    expectedStart = test_exp_df$Expected_Start,
    expectedRt = test_exp_df$Expected_Rt,
    expectedEnd = test_exp_df$Expected_End,
    verbose = TRUE
  )))
  
  expect_false(fail_result[[1]])
  
})
