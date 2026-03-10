library(LLSR)

test_that("search helpers return consistent CAS names", {
  cas_db <- LLSR::llsr_data$db.cas
  cas_search_db <- LLSR:::build_cas_search_db(cas_db)
  # Pick a known substring from the dataset (first CAS.NAME)
  first_name <- cas_db$CAS.NAME[1]
  substr_query <- substr(first_name, 1, min(5, nchar(first_name)))
  matches <- LLSR:::search_component_names(substr_query, cas_search_db)
  expect_true(first_name %in% matches)
})

test_that("AQSysCurve uses deterministic x-grid", {
  pars <- data.frame(42, -0.2, 1e-05)
  res1 <- AQSysCurve("merchuk", pars, silent = TRUE)
  res2 <- AQSysCurve("merchuk", pars, silent = TRUE)
  expect_equal(res1$data$X, res2$data$X)
})

