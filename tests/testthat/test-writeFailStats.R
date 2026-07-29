
make_fake_failure_df <- function(fake_dxf_folder, file_name, analysis, check_name) {
  file.create(file.path(fake_dxf_folder, file_name))
  data.frame(
    Analysis = analysis,
    fileID = file_name,
    failed_check = check_name,
    stringsAsFactors = FALSE
  )
}

make_fake_processed_list <- function(fake_dxf_folder, include_ref_samp_int = TRUE) {
  processed.list <- list(
    make_fake_failure_df(fake_dxf_folder, "failed_peak_number.dxf", 101, "peak_number"),
    list(make_fake_failure_df(fake_dxf_folder, "failed_reference_times.dxf", 102, "reference_times")),
    make_fake_failure_df(fake_dxf_folder, "failed_reference_isotope_sd.dxf", 103, "reference_isotope_sd"),
    make_fake_failure_df(fake_dxf_folder, "failed_reference_intensity.dxf", 104, "reference_intensity"),
    make_fake_failure_df(fake_dxf_folder, "failed_sample_processing.dxf", 105, "sample_processing"),
    make_fake_failure_df(fake_dxf_folder, "failed_sample_isotope_sd.dxf", 106, "sample_isotope_sd"),
    list(data.frame(Analysis = 201, fileID = "passed_refs.dxf")),
    list(data.frame(Analysis = 202, fileID = "passed_samps.dxf"))
  )

  if (include_ref_samp_int) {
    processed.list[[9]] <- make_fake_failure_df(
      fake_dxf_folder,
      "failed_reference_sample_intensity.dxf",
      107,
      "reference_sample_intensity"
    )
  }

  processed.list
}

test_that("writeFailStats writes every failure category", {
  oldwd <- getwd()
  test_dir <- tempfile("writeFailStats_all_failures")
  fake_dxf_folder <- file.path(test_dir, "fake_dxf_folder")
  output_dir <- file.path(test_dir, "csv_output")
  dir.create(fake_dxf_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(output_dir)
  on.exit(setwd(oldwd), add = TRUE)

  processed.list <- make_fake_processed_list(fake_dxf_folder, include_ref_samp_int = TRUE)

  result <- QCIRMS:::writeFailStats(processed.list = processed.list, outputID = "test")
  expected_failures <- data.frame(
    csv_file = c(
      file.path(output_dir, "test_failedPkNr.csv"),
      file.path(output_dir, "test_failedRefT.csv"),
      file.path(output_dir, "test_failedRefIso.csv"),
      file.path(output_dir, "test_failedRefInt.csv"),
      file.path(output_dir, "test_failedSampProc.csv"),
      file.path(output_dir, "test_failedSampIso.csv"),
      file.path(output_dir, "test_failedRefSampInt.csv")
    ),
    fileID = c(
      "failed_peak_number.dxf",
      "failed_reference_times.dxf",
      "failed_reference_isotope_sd.dxf",
      "failed_reference_intensity.dxf",
      "failed_sample_processing.dxf",
      "failed_sample_isotope_sd.dxf",
      "failed_reference_sample_intensity.dxf"
    ),
    failed_check = c(
      "peak_number",
      "reference_times",
      "reference_isotope_sd",
      "reference_intensity",
      "sample_processing",
      "sample_isotope_sd",
      "reference_sample_intensity"
    ),
    stringsAsFactors = FALSE
  )

  expect_null(result)
  expect_true(all(file.exists(expected_failures$csv_file)))
  expect_equal(length(list.files(fake_dxf_folder, pattern = "\\.dxf$")), 7)

  for (i in seq_len(nrow(expected_failures))) {
    csv_data <- read.csv(expected_failures$csv_file[i])

    expect_equal(csv_data$fileID, expected_failures$fileID[i])
    expect_true(file.exists(file.path(fake_dxf_folder, csv_data$fileID)))
    expect_equal(csv_data$failed_check, expected_failures$failed_check[i])
  }
})

test_that("writeFailStats skips reference-sample intensity output when absent", {
  oldwd <- getwd()
  test_dir <- tempfile("writeFailStats_no_ref_samp_int")
  fake_dxf_folder <- file.path(test_dir, "fake_dxf_folder")
  output_dir <- file.path(test_dir, "csv_output")
  dir.create(fake_dxf_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(output_dir)
  on.exit(setwd(oldwd), add = TRUE)

  processed.list <- make_fake_processed_list(fake_dxf_folder, include_ref_samp_int = FALSE)

  result <- QCIRMS:::writeFailStats(processed.list = processed.list, outputID = "test")
  expect_null(result)
  expect_false(file.exists(file.path(output_dir, "test_failedRefSampInt.csv")))
  expect_true(file.exists(file.path(output_dir, "test_failedSampIso.csv")))
})
