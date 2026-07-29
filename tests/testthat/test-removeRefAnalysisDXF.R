
test_that("function returns expected data type", {
  filtered.list <- vector("list", 8)
  filtered.list[[5]] <- data.frame(failed_sampProc_Analysis = 99)
  filtered.list[[6]] <- data.frame(failed_sampIso_Analysis = 99)
  filtered.list[[7]] <- list(data.frame(Analysis = 1, PeakNr = 1))
  filtered.list[[8]] <- list(data.frame(Analysis = 1, PeakNr = 2))

  result <- QCIRMS:::removeRefAnalysisDXF(filtered.list = filtered.list)

  expect_type(result, "list")
})
