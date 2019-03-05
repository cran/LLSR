####################################################################################################################
#' @import svDialogs
####################################################################################################################
#' @rdname AQSysCurve
#' @title This functions plot a curve based in the chosen model and its parameters.
#' @description The function returns a plot after using the parameters and model given by the user.
#' @details The function owns predefined set of equations that can be seen below and must be used, with adequated parameters,
#' to return a plot which represent the chosen model.
#' @export AQSysCurve
#' @param modelName Equation to be used: merchuk, murugesan [type:string]
#' @param modelPars Model's parameters [type::data.frame]
#' @param xlbl Plot's Horizontal axis label.
#' @param ylbl Plot's Vertical axis label.
#' 
#' @param col Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cex Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexlab Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexaxis Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexmain Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexsub Legacy from plot package. For more details, see \code{\link{plot.default}}
#' 
#' @param xmax Maximum value for the Horizontal axis' value (bottom-rich component)  [type:double]
#' @param NP Number of points used to build the fitted curve. Default is 100. [type:Integer]
#' @param save Save the generated plot in the disk using path and filename provided by the user. [type:Logical]
#' @param HR Adjust Plot's text to be compatible with High Resolution size [type:Logical]
#' @param seriesNames A list of sequential names which will identify each system provided by the user in the dataSET variable. [type:List]
#' @param filename Filename provided by the user to save a given plot. [type:String]
#' @param wdir The directory in which the plot file will be saved. [type:String]
#' @param silent save plot file without actually showing it to the user. [type:Logical]
# ' @param maxiter	- A positive integer specifying the maximum number of iterations allowed.
#' @inheritParams graphics::plot.default
#' @return A plot using the input model within the chosen interval and the curve's raw XY data.
#' If no interval is selected, xmax = 0.4.
#' @examples
#' \dontrun{
#' AQSysCurve("murugesan", data.frame(90.389, -34.897, 2.924), col = "red")
#' }
####################################################################################################################
AQSysCurve <- function  (
  modelName,
  modelPars,
  seriesNames = NULL,
  xlbl = "",
  ylbl = "",
  col = "black",
  type = "p",
  cex = 1,
  cexlab = 1,
  cexaxis = 1,
  cexmain = 1,
  cexsub = 1,
  xmax = 35,
  HR = FALSE,
  NP = 100,
  filename = NULL,
  wdir = NULL,
  save = FALSE,
  silent = FALSE,
  ...
)
{
  #
  nSys <- nrow(modelPars)
  models_npars <- AQSysList(TRUE)
  #
  if ((ncol(modelPars) == models_npars[[modelName]]) && (nSys >= 1 )) {
    # mass fraction range of bottom-rich component (min is 0, max is 1)
    x <- sort(runif(NP, 0.1, xmax))
    # select which model will be used to generate the plot
    # if user selects an option not available, it triggers an error
    # (check AQSys.err.R for details)
    Fn <- ifelse(modelName %in% names(models_npars), AQSys.mathDesc(modelName), AQSys.err("0"))
    #
    if (is.null(seriesNames) || !(length(seriesNames) == nSys)) {
      print(paste("The array seriesNames must have", nSys, "elements. Default names will be used instead."))
      seriesNames <- sapply(seq(1, nSys), function(x) paste("Series", x))
    } else {
      SysNames <- TRUE
    }
    # CREATE A LIST WHICH WILL HOLD EACH SYSTEM'S DATA
    SysList <- list()
    for (i in seq(1, nSys)) {
      # unlist and convert parameters to double
      model_pars <- as.double(unlist(modelPars[i, ]))
      SysList[[i]] <- unname(data.frame(x, Fn(model_pars, x)))
      names(SysList[[i]]) <- c("X", "Y")
      SysList[[i]]["System"] <- seriesNames[i]
    }
    # BIND ALL DATA FROM SEVERAL SYSTEMS
    output_data <- bind_rows(SysList)
    # PLOT DATA
    output_plot <- BLOPlot(output_data, xlbl, ylbl)
    #
    saveConfig(output_plot, save, HR, filename, wdir, silent)
    # make available data from fitted curve to user. Function returns it silently
    # but user can get data using simple assign '<-'
    if (silent == FALSE) {
      print(output_plot)
      invisible(output_data)
    } else {
      invisible(list("data" = output_data, "plot" = output_plot))
    }
  }
  else{
    AQSys.err("9")
  }
}

BLOPlot <- function(dataSET, xlbl = "", ylbl = "") {
  # Calculate X and Y maximums as multiple of five and 8% superior of the axis maximum value
  xmax <- ceiling(round(max(dataSET$X) / 0.92, 1) / 5) * 5
  ymax <- ceiling(round(max(dataSET$Y) / 0.92, 1) / 5) * 5
  # Plot data
  output_plot <- ggplot(dataSET, aes_string(x = "X", y = "Y", color = "System", shape = "System"))  +
    geom_line(size = 1) + 
    geom_point(size = 2) +
    theme_light() +  
    scale_color_llsr(palette = "mixed") +
    xlab(paste(xlbl,  "(%, m/m)")) +
    ylab(paste(ylbl, "(%, m/m)")) + 
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
      limits = c(-2.5, ymax),
      breaks = seq(0, ymax, by = 5),
      labels = seq(0, ymax, by = 5)
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(-0.5, xmax),
      breaks = seq(0, xmax, by = xmax / 10),
      labels = seq(0, xmax, by = xmax / 10)
    )
  return(output_plot)
}
