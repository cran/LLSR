library(LLSR)
options(digits = 15)

test_that("Mechuk's Coefficients", {
  expect_equal(summary(AQSys(llsr_data$db.data[6:23, 1:2], 
                             modelName = "merchuk"))$coefficients[1], 
               42.9217635287213, tolerance = 1e-7)
  expect_equal(summary(AQSys(llsr_data$db.data[6:23, 1:2], 
                             modelName = "merchuk"))$coefficients[2],
               -0.229660343210977,
               tolerance = 1e-7)
  expect_equal(summary(AQSys(llsr_data$db.data[6:23, 1:2], 
                             modelName = "merchuk"))$coefficients[3],
               1.45845232545921e-05,
               tolerance = 1e-7)
})

test_that("Murugesan's Coefficients", {
  expect_equal(summary(AQSys(llsr_data$db.data[6:23, 1:2], 
                             modelName = "murugesan"))$coefficients[1], 
               37.9239775477034, tolerance =
                 1e-7)
  expect_equal(summary(AQSys(llsr_data$db.data[6:23, 1:2], 
                             modelName = "murugesan"))$coefficients[2], 
               -5.64719955175667, tolerance =
                 1e-7)
  expect_equal(summary(AQSys(llsr_data$db.data[6:23, 1:2],
                             modelName = "murugesan"))$coefficients[3],
               0.0450678037108678,
               tolerance = 1e-7)
})

# AQSys.gsnchk(peg4kslt[1:2],peg4kslt[2,3],peg4kslt[2,4],peg4kslt[2,5],
# peg4kslt[2,6],peg4kslt[2,7],peg4kslt[2,8])
