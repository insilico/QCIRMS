
# make_brines_fixture <- function(n = 1) {
#   brines_dir <- normalizePath(
#     file.path(testthat::test_path(), "..", "..", "220318_1pct_CO2_NaCl_NaSO4_brines"),
#     mustWork = TRUE
#   )
#   dxf_files <- list.files(brines_dir, pattern = "\\.dxf$", full.names = TRUE)
#   fixture_dir <- file.path(tempdir(), paste0("qcirms_brines_", Sys.getpid(), "_", sample.int(1e6, 1)))
#   dir.create(fixture_dir)
#   file.copy(dxf_files[seq_len(n)], fixture_dir, overwrite = TRUE)
#   fixture_dir
# }
# 

test_that("function returns expected data type", {

  result <- find_files_by_analysis_num(path = builtin_dxf_folder, analysisNum.vec = test_analysisNum.vec)
  
  expect_type(result, "character")
})

