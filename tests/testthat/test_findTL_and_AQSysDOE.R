library(LLSR)

test_that("findTL errors clearly when requested TLL is out of bounds", {
  # Build a minimal synthetic SysTLL structure consistent with expected fields
  minTL <- data.frame(X = c(1, 2, 3), Y = c(1, 2, 3))
  maxTL <- data.frame(X = c(4, 5, 6), Y = c(4, 5, 6))
  SysTLL <- list(
    TLL = TLL(minTL, maxTL),
    maxTL = maxTL
  )
  BLFn <- function(x) x  # simple diagonal
  # dTLL greater than MaxTLL should trigger an informative error
  expect_error(
    findTL(SysTLL$TLL$MaxTLL * 1.5, SysTLL, BLFn, slope = 1),
    "Target-TLL must be BIGGER"
  )
})


