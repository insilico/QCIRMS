
test_that("function returns expected data type", {
  allStandards_d18O.list <- list(
    data.frame(Analysis = 1, Identifier1 = "L1", PkNr = 1:3, d18O16O = c(1, 1.1, 1.2), d13C12C = 1),
    data.frame(Analysis = 2, Identifier1 = "H1", PkNr = 1:3, d18O16O = c(2, 2.1, 2.2), d13C12C = 1),
    data.frame(Analysis = 3, Identifier1 = "LW", PkNr = 1:3, d18O16O = c(3, 3.1, 3.2), d13C12C = 1)
  )

  result <- QCIRMS:::avg_sd_d18O_standards(
    allStandards_d18O.list = allStandards_d18O.list,
    standNames = c("L1", "H1", "LW"),
    standAcceptedVals_vec = c(1, 2, 3),
    accStandRatioSD = c(0.2, 0.2, 0.2)
  )

  expect_type(result, "list")
})
