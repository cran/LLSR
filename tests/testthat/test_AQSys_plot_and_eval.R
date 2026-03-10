library(LLSR)

test_that("AQSys.plot uses xmax for x-axis scale", {
  dataSET <- llsr_data$db.data[6:23, 1:2]
  p <- AQSys.plot(dataSET, xmax = 50, ymax = 30, silent = TRUE)
  # AQSys.plot returns the plot object when silent = TRUE
  expect_s3_class(p, "ggplot")
  panel <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
  # Be robust across ggplot2 versions: prefer x.range when available,
  # otherwise fall back to x$range or the range of breaks.
  x_limits <- if (!is.null(panel$x.range)) {
    panel$x.range
  } else if (!is.null(panel$x$range)) {
    panel$x$range
  } else {
    range(panel$x$breaks)
  }
  expect_equal(x_limits, c(0, 50))
})

test_that("AQSysEval returns matching data and plot structures", {
  dataSET <- llsr_data$db.data[6:23, 1:2]
  # Provide a dummy slope to avoid dependence on db.tielines contents in tests
  eval_res <- AQSysEval(
    dataSET,
    xlbl = names(dataSET)[1],
    ylbl = names(dataSET)[2],
    slope = 1
  )
  # Single-system call should return a list for 'data' and a ggplot for 'plot'
  expect_true(is.list(eval_res$data))
  expect_s3_class(eval_res$plot, "ggplot")
})

