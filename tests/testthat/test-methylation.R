context("Processing data")

test_that("data centering works", {
  region <- c(2, 4, 5)
  tss <- 3
  expect_identical(center_loc(region, tss, "+"), c(-1, 1, 2))
  expect_identical(center_loc(region, tss, "-"), c(1, -1, -2))
  expect_error(center_loc(region, tss, 9))
  expect_error(center_loc(region, "tss", 9))
  expect_identical(center_loc(21, 2, "+"), 19)
})

test_that("minmax scaling works",{
  data <- c(-20, 0, 15, 20)
  xmin <- -20
  xmax <- 20
  expect_identical(minmax_scaling(data, xmin, xmax), c(0, 0.5, 0.875, 1))
  expect_error(minmax_scaling(data, -xmin, xmax))
})
