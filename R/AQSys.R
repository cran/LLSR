if(getRversion() >= "3.5")
  utils::globalVariables(
    c(
      "x",
      "Series",
      "y",
      "xlab",
      "ylab",
      "geom_point",
      "geom_line",
      "unit",
      "element_text",
      "element_line",
      "element_rect",
      "scale_x_continuous",
      "scale_y_continuous"
    )
  )
####################################################################################################################
#' @import rootSolve
#' @import graphics
#' @import stats
#' @import svDialogs
#' @import ggplot2
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
####################################################################################################################
require(rootSolve)
require(svglite)
require(svDialogs)
require(ggplot2)
# options(warn = 1)
####################################################################################################################
#'Merchuk's Equation to fit Binodal Experimental Data
#' @rdname AQSys
#' @name AQSys
#' @description .
#' @export
#' @seealso \itemize{
#' \item \code{\link{AQSys.default}}
#' \item \code{\link{AQSys.plot}}
#' \item \code{\link{AQSys.tielines}}
#' \item \code{\link{AQSys.crpt}}
#' \item \code{\link{AQSysOthmer}}
#' \item \code{\link{AQSysBancroft}}
#' }
####################################################################################################################
AQSys <- function(XYdt, ...)
  UseMethod("AQSys")
####################################################################################################################
#' @rdname AQSys
#' @title Merchuk's nonlinear Equation
#' @description Perform a nonlinear regression fit using several mathemmatical descriptors in order to determine the
#' equation's parameters.
#' @details The function returns functions parameters after fitting experimental data to the equations listed in AQSysList().
#' @param mathDesc - Character String specifying the nonlinear empirical equation to fit data.
#' The default method uses Merchuk's equation. Other mathematical descriptors can be listed using AQSysList().
#' @param XYdt - Binodal Experimental data that will be used in the nonlinear fit
#' @param ... Additional optional arguments. None are used at present.
#' @method AQSys default
#' @export
#' @return A list containing three data.frame variables with all data parsed from the worksheet and parameters calculated
#' through the available mathematical descriptions.
#' @examples
#' #Populating variable XYdt with binodal data
#' XYdt <- peg4kslt[,1:2]
#' #Fitting XYdt using Merchuk's function
#' AQSys(XYdt)
AQSys.default <- function(XYdt, mathDesc = "merchuk", ...) {
  # each switch option calls a correspondent equation to fit XYdt
  # equations are functions declared in AQSysFormulas.R
  switch(
    mathDesc,
    merchuk = {
      ans <- mrchk(XYdt)
    },
    murugesan = {
      ans <- mrgsn(XYdt)
    },
    tello = {
      ans <- tello(XYdt)
    },
    AQSys.err("0")
  )
  # return fitting parameters ans statistical data
  return(ans)
}
#
####################################################################################################################
# MERCHUK PLOT TEST FUNCTION
#' @rdname AQSys.plot
#' @title Dataset and Fitted Function plot
#' @description The function returns a plot after fitting a dataset to a given equation.
#' @details This version uses the plot function and return a regular bidimensional plot.
#' @method AQSys plot
#' @export AQSys.plot
#' @export
#' @param mathDesc - Character String specifying the nonlinear empirical equation to fit data. The default method uses
#' Merchuk's equation. Other possibilities can be seen in AQSysList().
#' @param ... Additional optional arguments. None are used at present.
#' @param XYdt Binodal Experimental data that will be used in the nonlinear fit
#' @param xlbl Plot's Horizontal axis label.
#' @param ylbl Plot's Vertical axis label.
#' @param main Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param col Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param type Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cex Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexlab Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexaxis Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexmain Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexsub Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param xmax Maximum value for the Horizontal axis' value
#' @param ymax Maximum value for the Vertical axis' value
#' @param filename A filename chosen by the user to save a given plot
#' @param wdir A directory in which the plot file will be saved
#' @param silent save the file without the user input
#' @param HR Magnify Plot's text to be compatible with High Resolution size [type:Boulean]
#' @param NP Number of points used to build the fitted curve. Default is 100. [type:Integer]
#' @param clwd Plot's axis line width
#' @param save Optimize Plot's elements to be compatible with High Resolution size [type:Boulean]
#' @return A plot containing the experimental data, the correspondent curve for the binodal in study and the curve's raw XY data.
#' @examples
#' #Populating variable XYdt with binodal data
#' XYdt <- peg4kslt[,1:2]
#' #Plot XYdt using Merchuk's function
#' #
#' AQSys.plot(XYdt)
#' #
AQSys.plot <-
  function  (XYdt,
             xlbl = "",
             ylbl = "",
             main = NULL,
             col = "blue",
             type = "p",
             cex = 1,
             cexlab = 1,
             cexaxis = 1,
             cexmain = 1,
             cexsub = 1,
             xmax = "",
             ymax = "",
             HR = FALSE,
             NP = 100,
             mathDesc = "merchuk",
             clwd = NULL,
             filename = NULL,
             wdir = NULL,
             save = FALSE,
             silent = FALSE,
             ...)
  {
    #
    plot_image = NULL
    #
    if (xlbl == "") {
      xlbl <- names(XYdt)[1]
    }
    if (ylbl == "") {
      ylbl <- names(XYdt)[2]
    }
    names(XYdt) <- c('x', 'y')
    #
    if ((xmax == "") ||
        (xmax > 1) || (xmax < round(max(XYdt[1]) / 0.92, 1))) {
      xmax <- round(max(XYdt[1]) / 0.92, 1)
    }
    if ((ymax == "") ||
        (ymax > 1) || (ymax < round(max(XYdt[2]) / 0.92, 1))) {
      ymax <- round(max(XYdt[2]) / 0.92, 1)
    }
    #
    #
    if (save == TRUE) {
      if (HR == TRUE) {
        image_format <- ".svg"
      } else{
        image_format <- ".png"
      }
      #
      if (is.null(filename)) {
        # Get user choice for a filename to save the plot
        filename <-
          dlgInput(message = "Enter the figure filename:")$res
      }
      # complete filename with the appropriated extension
      filename <- paste(filename, image_format, sep = "")
      # Check if filename is invalid and quite if so
      if (filename == image_format) {
        stop("Filename is NULL or INVALID.", call. = TRUE)
      }
      #
      #
      if (is.null(wdir)) {
        # Get user choice for a directory to save the plot
        wdir <- dlgDir()$res
      }
      # Check if path is invalid and quite if so
      if ((wdir == "") && (silent == FALSE)) {
        #
        stop("Path is NULL or INVALID.", call. = TRUE)
        #
      } else if ((wdir == "") && (silent == TRUE)) {
        #
        wdir <- getwd()
        wdir <- paste(wdir, filename, sep = .Platform$file.sep)
        #
      } else{
        #
        wdir <- paste(wdir, filename, sep = .Platform$file.sep)
        #
      }
    }
    #
    #
    #
    # Select which model will be used to generate the plot
    # SWITCH WORKS ONLY FOR THREE PARAMETER'S EQUATIONS.
    # IF NECESSARY A HIGHER NUMBER, INSERT CONDITIONAL BELOW.
    switch(
      mathDesc,
      merchuk = {
        # fit data using chosen equation and get coefficients
        CoefSET <- summary(mrchk(XYdt))$coefficients[, 1]
        # set the equation that will be used to plot the phase diagram
        Fn <- AQSys.mathDesc("merchuk")
      },
      murugesan = {
        # fit data using chosen equation and get coefficients
        CoefSET <- summary(mrgsn(XYdt))$coefficients[, 1]
        # set the equation that will be used to plot the phase diagram
        Fn <- AQSys.mathDesc("murugesan")
      },
      tello = {
        # fit data using chosen equation and get coefficients
        CoefSET <- summary(tello(XYdt))$coefficients[, 1]
        # set the equation that will be used to plot the phase diagram
        Fn <- AQSys.mathDesc("tello")
      },
      # if user selects an option not available, it triggers an error (check AQSys.err.R for details)
      AQSys.err("0")
    )
    #
    #
    #
    #plot phase diagram using experimental data and with previously calculated parameters
    plot_image <-
      ggplot(data = XYdt, aes(x, y)) + geom_point(shape = 8, size = 2) +
      theme_light() + xlab(xlbl) + ylab(ylbl) + theme(validate = FALSE,
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(size = 18),
        legend.position = "top",
        axis.title.y = element_text(vjust = 1.4),
        axis.title.x = element_text(vjust = 0),
        panel.grid.major = element_line(size = .70, colour = "black"),
        panel.grid.minor = element_line(size = .70),
        panel.border = element_rect(size = .5, colour = "white"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.line = element_line(colour = 'black', size = 1.25),
        legend.title = element_text(
          colour = "black",
          size = 12,
          face = "bold",
          angle = 0
        ),
        legend.text = element_text(
          colour = "black",
          size = 12,
          face = "plain"
        )
      ) +
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0.001, ymax),
        breaks = seq(0, ymax, by = ymax / 10),
        labels = seq(0, ymax, by = ymax / 10)
      ) +
      scale_x_continuous(
        expand = c(0, 0),
        limits = c(0, xmax),
        breaks = seq(0, ymax, by = xmax / 10),
        labels = seq(0, ymax, by = xmax / 10)
      )
    #
    # add curve generated using regression parameters
    rawdt <- data.frame(XYdt[1], Fn(CoefSET, XYdt[1]))
    names(rawdt) <- c("x", "y")
    plot_image <-
      plot_image + geom_line(
        data = rawdt,
        aes(x = x, y = y),
        color = "red",
        linetype = 2,
        size = 1.1
      )
    #
    if (save == TRUE) {
      ggsave(
        filename = wdir,
        plot = plot_image,
        width = 21.14 / 2,
        height = 14.39 / 2
      )
    }
    #
    # make available data from fitted curve to user. Function returns it silently
    # but user can get data using simple assign '<-'
    #
    if (silent == FALSE) {
      print(plot_image)
      invisible(rawdt)
    } else {
      invisible(plot_image)
    }
    #
  }
####################################################################################################################
