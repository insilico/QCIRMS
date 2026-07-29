
test_that("function returns expected data type", {
  result <- peak_area_trap(
    start.t = 1,
    end.t = 3,
    time.vec = 1:6,
    int.vec = c(0, 2, 0, 0, 3, 0)
  )

  expect_type(result, "double")
})
