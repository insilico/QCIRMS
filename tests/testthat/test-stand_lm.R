
test_that("function returns expected data type", {
  acceptedMeas.df <- data.frame(
    accepted_d18O16O = c(1, 2, 3),
    measured_d18O16O = c(1.1, 2, 3.4),
    accepted = c(1, 2, 3),
    measured = c(1.1, 2, 3.4)
  )
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  result <- QCIRMS:::stand_lm(acceptedMeas.df = acceptedMeas.df, dataName = "test")

  expect_type(result, "list")
})
