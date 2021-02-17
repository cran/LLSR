###############################################################################
options(digits = 14)
###############################################################################
#' @import ggplot2
#' @import rootSolve
###############################################################################
crit_point_poly <- function(dataSET,
                            tldt,
                            modelName,
                            xmax,
                            xlbl,
                            ylbl,
                            Order,
                            ext) {
  #
  required_fields <- c("A", "B", "ORDER", "PH", "TEMP", "TOP.A", "TOP.B", 
                       "BOT.A", "BOT.B", "GLB.A", "GLB.B")
  #
  poly_data <- setNames(data.frame(matrix(nrow = 0, ncol = 2)), 
                        c("X", "TLL"))
  dataSET <- suppressWarnings(toNumeric(dataSET, Order))
  #
  if (all(required_fields %in% names(tldt))) {
    dataTL <- setNames(data.frame(matrix(nrow = 0, ncol = 3)),
                       c("X", "Y", "TL"))
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
      Ym <- sum(Ys) / 2
      dX <- diff(Xs)
      Xm <- sum(Xs) / 2
      #
      tempTL <- data.frame(X = c(Xs, Xm), Y = c(Ys, Ym), TL = rep(row, 3))
      dataTL <- rbind(dataTL, tempTL)
      #
      slope <- (dY / dX)
      #
      tll <- sqrt((dX ^ 2) + (dY ^ 2))
      #
      row_entry <- data.frame(X = Xm, TLL = tll)
      poly_data <- rbind(poly_data, row_entry)
    }
    #
    rownames(poly_data) <- NULL
    poly_model <- lm(poly_data$TLL ~ poly(poly_data$X, 3, raw = TRUE))
    # n <- summary(model)
    # n$r.squared
    coefs <- unname(unlist(lapply(poly_model$coefficients, function(ith_coeff) {
      ifelse(is.na(ith_coeff), 0, ith_coeff)
    })))
    #
    BNNLAnalysis <- AQSys(dataSET, modelName = modelName)
    PARs <- t(summary(BNNLAnalysis)$parameters[, 1])
    #
    xmax <- ifelse((xmax == "" | is.null(xmax)), 
                   ceiling(round(max(dataSET[, 1]) / 0.92, 1) / 5) * 5, xmax)
    #
    BNFNs <- mathDescPair(modelName)
    #
    EqSys <- function(x) {
      F1 <- eval(parse(text = gsub("[$]", "", BNFNs)))
      F2 <- coefs[1] + coefs[2]*x[2] +
        coefs[3]*(x[2]^2) + 
        coefs[4]*(x[2]^3)-x[1]
      #
      c(F1 = F1, F2 = F2)
    }
    OUTPUT <- setNames(multiroot(
      f = EqSys,
      start = c(10, 10),
      positive = TRUE
    )$root, c("YC", "XC"))
    #
    if (ext) {
      x <- sort(runif(100, 0.1, xmax))
      Y <- function(x) {
        coefs[1] + coefs[2] * x + coefs[3] * (x ^ 2) + coefs[4] * (x ^ 3)
      }
      xy <- data.frame(Xs = x, Ys = Y(x))
      #
      BNPlot <- AQSys.plot(dataSET,
                           silent = TRUE,
                           xmax = xmax,
                           xlbl = xlbl,
                           ylbl = ylbl)
      PolyPlot <- BNPlot + 
        # geom_line(
        #   data = xy,
        #   aes_string(x = "Xs", y = "Ys"),
        #   size = 1.1,
        #   linetype = "dashed",
        #   color = "cornflowerblue"
        # ) +
        geom_line(
          data = dataTL,
          aes_string(x = "X", y = "Y", group = "TL"),
          colour = "red",
          alpha = 0.4
        ) +
        geom_point(
          data = dataTL,
          aes_string(x = "X", y = "Y", group = "TL"),
          colour = "black",
          bg = "red",
          shape = 21,
          alpha = 1
        )
      OUTPUT_PLOT <- PolyPlot + annotate(
        "point",
        x = OUTPUT[2],
        y = OUTPUT[1],
        colour = "black",
        bg = "gold",
        shape = 23,
        size = 2
      )
      # return(dataTL)
      return(list(CriticalPoint=OUTPUT[2:1], Plot=OUTPUT_PLOT))
    }
    return(OUTPUT[2:1])
  } else {
    # trigger error
  }
}
###############################################################################
