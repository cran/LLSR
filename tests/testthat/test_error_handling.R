library(LLSR)

test_that("AQSys.default validates input type and shape", {
  # Non-data.frame input should error
  expect_error(AQSys(matrix(1:4, ncol = 2)), "must be a data.frame")
  # Data.frame with fewer than 2 columns should error
  df_bad <- data.frame(x = 1:3)
  expect_error(AQSys(df_bad), "valid data set has multiple two-columns data")
})

test_that("AQSearch requires at least one search parameter", {
  expect_error(
    AQSearch(),
    "At least one of the parameters"
  )
})

