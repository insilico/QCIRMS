
test_that("function reports the current generic_plot_all_raw argument mismatch", {
  expect_error(
    plot_pass_fail_spectra(
      samps.dat = filtered_check_df,
      dataName = "test",
      plotsPath = tempdir(),
      dxf_path = builtin_dxf_folder
    ),
    "unused argument"
  )
})
