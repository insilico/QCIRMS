
test_that("function reports the current generic_plot_all_raw argument mismatch", {

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_error(
    plot_by_analysis_num(
      dxf_path = builtin_dxf_folder,
      analysisNum.vec = test_analysisNum.vec,
      write_pdf = FALSE,
      plot_path = tempdir(),
      pdf_name = "test_spectra.pdf"
    )
  )
})
