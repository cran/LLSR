#' @importFrom dplyr bind_rows
####################################################################################################################
#' @rdname AQSysPlot
#' @title This functions plot a curve based in the chosen model and its parameters.
#' @description The function returns a plot after using the parameters and model given by the user.
#' @details The function owns predefined set of equations that can be seen below and must be used, with adequated parameters,
#' to return a plot which represent the chosen model.
#' 
#' @param dataSET - Binodal Experimental data that will be used in the nonlinear fit. It might hold multiple systems stacked side-by-side. [type:data.frame]
#' @param xlbl Plot's Horizontal axis label.
#' @param ylbl Plot's Vertical axis label.
#' @param seriesNames A list of sequential names which will identify each system provided by the user in the dataSET variable. [type:List]
#' @param save Save the generated plot in the disk using path and filename provided by the user. Default is False. [type:Logical]
#' @param filename A filename provided by the user to save a given plot. No default is provided. [type:String]
#' @param HR  Adjust Plot's text to be compatible with High Resolution size [type:Logical]
#' @param wdir The directory in which the plot file will be saved. [type:String]
#' @param silent save plot file without actually showing it to the user. Default is FALSE. [type:Logical]
#' @export AQSysPlot
AQSysPlot <- function (dataSET,
                       xlbl,
                       ylbl,
                       seriesNames = NULL,
                       save = FALSE,
                       filename = NULL,
                       HR = FALSE,
                       wdir = NULL,
                       silent = FALSE) {
  #
  nSys <- (ncol(dataSET) / 2)
  SysNames <- FALSE
  #
  if ((ncol(dataSET) %% 2) == 0) {
    if (is.null(seriesNames) || !(length(seriesNames) == nSys)) {
      print(paste("The array seriesNames must have", nSys, "elements. Default names will be used instead."))
      seriesNames <- sapply(seq(1, nSys), function(x) paste("Series", x))
    } else {
      SysNames <- TRUE
    }
    SysList <- list()
    for (i in seq(1, nSys)) {
      SysList[[i]] <- unname(na.exclude(dataSET[, (i * 2 - 1):(i * 2)]))
      names(SysList[[i]]) <- c("X", "Y")
      SysList[[i]]["System"] <- seriesNames[i]
    }
    output <- bind_rows(SysList)
    output_plot <- bndOrthPlot(output, xlbl, ylbl)
    #
    saveConfig(output_plot, save, HR, filename, wdir, silent)
    #
    if (silent == FALSE) {
      print(output_plot)
      invisible(list("data" = SysList, "plot" = output_plot))
    } else {
      invisible(list("data" = SysList, "plot" = output_plot))
    }
  }
  else{
    AQSys.err(9)
  }
}

bndOrthPlot <- function(dataSET, xlbl = "", ylbl = "") {
  #
  xmax <- ceiling(round(max(dataSET[, 1]) / 0.92, 1) / 5) * 5
  ymax <- ceiling(round(max(dataSET[, 2]) / 0.92, 1) / 5) * 5
  #
  outputPLOT <- ggplot() + scale_colour_grey() +
    geom_line(data=dataSET, size = 1, aes_string(color = "System", x = "X", y = "Y")) +
    geom_point(data=dataSET, size = 2, aes_string(color = "System", x = "X", y = "Y")) + 
    xlab(paste(xlbl,  "(%, m/m)")) +
    ylab(paste(ylbl, "(%, m/m)")) + 
    theme_light() + 
    theme(
      validate = FALSE,
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      text = element_text(size = 16),
      legend.position = "top",
      axis.title.y = element_text(vjust = 5),
      axis.title.x = element_text(vjust = -2),
      panel.grid.major = element_line(size = .70, colour = "black"),
      panel.grid.minor = element_line(size = .70),
      panel.border = element_rect(size = .5, colour = "white"),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      axis.line = element_line(colour = 'black', size = 1.25),
      legend.title = element_blank(),
      legend.text = element_text(
        colour = "black",
        size = 12,
        face = "plain"
      )
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(-5, ymax),
      breaks = seq(0, ymax, by = 5),
      labels = seq(0, ymax, by = 5)
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, xmax),
      breaks = seq(0, xmax, by = xmax / 10),
      labels = seq(0, xmax, by = xmax / 10)
    )
  #
  return(outputPLOT)
}


