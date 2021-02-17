###############################################################################
options(digits = 14)
###############################################################################
#' @import ggplot2
#' @import nleqslv
#' @import rootSolve
###############################################################################
crit_point_eqsys <- function(dataSET,
                             tldt,
                             modelName,
                             xmax,
                             xlbl,
                             ylbl,
                             Order,
                             ext) {
  
  ###
  dataSET <- suppressWarnings(toNumeric(dataSET, Order))
  PARs <- summary(AQSys(dataSET))$coefficients[, 1]
  BnFn <- mathDescPair(modelName)
  d <- 10^-100
  #
  #
  #
  poly_data <- setNames(data.frame(matrix(nrow = 0, ncol = 2)), c("S", "TLL"))
  #
  for (row in seq(1, nrow(tldt))) {
    tldt_row <- tldt[row, ]
    #
    if (tolower(tldt_row["ORDER"]) == "yx") {
      Ys <- unlist(tldt_row[c("TOP.A", "BOT.A")])
      Xs <- unlist(tldt_row[c("TOP.B", "BOT.B")])
    } else {
      Xs <- unlist(tldt_row[c("TOP.A", "BOT.A")])
      Ys <- unlist(tldt_row[c("TOP.B", "BOT.B")])
    }
    #
    dY <- diff(Ys)
    dX <- diff(Xs)
    slope <- (dY / dX)
    tll <- sqrt((dX ^ 2) + (dY ^ 2))
    #
    row_entry <- data.frame(S = slope, TLL = tll)
    poly_data <- rbind(poly_data, row_entry)
  }
  rownames(poly_data) <- NULL
  poly_model <- lm(poly_data$S ~ poly(poly_data$TLL, 3, raw = TRUE))
  coefs <- unname(unlist(lapply(poly_model$coefficients, function(ith_coeff) {
    ifelse(is.na(ith_coeff), 0, ith_coeff)
  })))
  SFn <- function(x) {
    sum(coefs[1] +
          coefs[2] * x +
          coefs[3] * (x ^ 2) +
          coefs[4] * (x ^ 3))
  }
  m <- SFn(0)
  #
  #
  #
  EqSys <- function(x) {
    F1 <- eval(parse(text = gsub("[$]", "", BnFn)))
    F2 <- m*(x[2] - x[4]) - (x[1] - x[3])
    F3 <- ((x[2] - x[4]) ^ 2) + ((x[1] - x[3]) ^ 2) - d^2
    F4 <- (((0 - x[4]) ^ 2) + ((0 - x[3]) ^ 2)) - 
      (((x[2] - 0) ^ 2) + ((x[1] - 0) ^ 2))
    #
    return(c(F1 = F1,
             F2 = F2,
             F3 = F3,
             F4 = F4))
  }
  #
  sysres <- nleqslv(c(10, 10, 10 ^ -10, 10), EqSys)$x
  OUTPUT <- setNames(sysres[1:2], c("YC", "XC"))
  #
  #
  #
  if (ext) {
    #
    OUTPUT_PLOT <- AQSys.plot(dataSET = dataSET,
                         silent = TRUE,
                         xmax = xmax,
                         xlbl = xlbl,
                         ylbl = ylbl) + annotate(
                           "point",
                           x = OUTPUT[2],
                           y = OUTPUT[1],
                           colour = "black",
                           bg = "gold",
                           shape = 23,
                           size = 2
                         )
    return(list(CriticalPoint=OUTPUT[2:1], Plot=OUTPUT_PLOT))
  }
  return(OUTPUT[2:1])
}
###############################################################################
