
test_that("function returns expected data type", {
  result <- all_PA_trap(
    start.vec = c(1, 4),
    end.vec = c(3, 6),
    time.vec = 1:6,
    int.vec = c(0, 2, 0, 0, 3, 0),
    pk.Nrs = c(1, 2)
  )

  expect_s3_class(result, "data.frame")
})
