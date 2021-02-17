if (getRversion() >= "3.5")
  utils::globalVariables(
    c(
      "Series",
      "wt",
      "wb",
      "xlab",
      "ylab",
      "geom_point",
      "geom_line",
      "unit",
      "element_text",
      "element_line",
      "element_rect",
      "scale_x_continuous",
      "scale_y_continuous",
      "A",
      "B",
      "PH",
      "TEMP"
    )
  )
###############################################################################
#' @import rootSolve graphics stats svDialogs ggplot2 svglite
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
# options(warn = 1)
###############################################################################
# Merchuk's Equation to fit Binodal Experimental Data
#' @rdname AQSys
#' @name AQSys
#' @description .
#' @export
#' @seealso \itemize{
#' \item \code{\link{AQSys.default}}
#' \item \code{\link{AQSys.plot}}
#' \item \code{\link{AQSys.LevArmRule}}
#' \item \code{\link{AQSysOthmer}}
#' \item \code{\link{AQSysBancroft}}
#' }
###############################################################################
AQSys <- function(dataSET, ...)
  UseMethod("AQSys")
###############################################################################
#' @rdname AQSys
#' @title Merchuk's nonlinear Equation
#' @description Perform a nonlinear regression fit using any of the several 
#' mathematical descriptors implemented in order to calculate the equation's 
#' parameters.
#' @details The function returns functions parameters after fitting experimental
#'  data to the equations listed in AQSysList().
#' @param modelName - Character String specifying the nonlinear empirical 
#' equation to fit data.
#' The default method uses Merchuk's equation. Other mathematical descriptors 
#' can be listed using AQSysList().
#' @param dataSET - Binodal Experimental data that will be used in the 
#' nonlinear fit
#' @param Order Defines how the data is organized in the Worksheet. Use "xy"
#'  whether the first column corresponds to the lower phase fraction and "yx"
#'   whether the opposite.
# ' @param maxiter	- A positive integer specifying the maximum number of 
# iterations allowed.
#' @param ... Additional optional arguments. None are used at present.
#' @method AQSys default
#' @export
#' @return A list containing three data.frame variables with all data parsed 
#' from the worksheet and parameters calculated
#' through the available mathematical descriptions.
#' @examples
#' # Populating variable dataSET with binodal data
#' dataSET <- peg4kslt[ , 1:2]
#' # Fitting dataSET using Merchuk's function
#' AQSys(dataSET)
#' @references 
#' MURUGESAN, T.; PERUMALSAMY, M. Liquid-Liquid Equilibria of 
#' Poly(ethylene glycol) 2000 + Sodium Citrate + Water at 
#' (25, 30, 35, 40, and 45) C. Journal of Chemical & Engineering Data,
#'  v. 50, n. 4, p. 1392-1395, 2005/07/01 2005. ISSN 0021-9568. 
#' (\doi{10.1021/je050081k})
#' 
#' MERCHUK, J. C.; ANDREWS, B. A.; ASENJO, J. A. Aqueous two-phase systems for 
#' protein separation: Studies on phase inversion. Journal of Chromatography B: 
#' Biomedical Sciences and Applications, v. 711, n. 1-2, p. 285-293,  1998. 
#' ISSN 0378-4347. (\doi{10.1016/s0378-4347(97)00594-x})
#' 
#' TANG, X.  et al. The study of phase behavior of aqueous two-phase system 
#' containing [Cnmim] BF 4 (n=2, 3, 4)+(NH4)2SO4 + H2O at different 
#' temperatures. Fluid Phase Equilibria, v. 383, p. 100-107,  2014. 
#' ISSN 0378-3812. (\doi{10.1016/j.fluid.2014.09.029})
#' 
#' GONZALEZ-TELLO, P.  et al. Liquid-Liquid Equilibrium in the System 
#' Poly(ethylene glycol) + MgSO4 + H2O at 298 K. Journal of Chemical &
#'  Engineering Data, v. 41, n. 6, p. 1333-1336, 1996/01/01 1996. 
#'  ISSN 0021-9568. (\doi{10.1021/je960075b})
#' 
#' CHEN, Y.  et al. Liquid-liquid equilibria of aqueous biphasic systems 
#' composed of 1-butyl-3-methyl imidazolium tetrafluoroborate+ sucrose/maltose+
#' water. Journal of Chemical & Engineering Data, v. 55, n. 9, p. 3612-3616,
#' 2010. ISSN 0021-9568. (\doi{10.1021/je100212p})
# 
# XUEQIAO, X.  et al. Measurement and Correlation of the Phase Diagram Data for
# PPG400 + (K3PO4, K2CO3, and K2HPO4) + H2O Aqueous Two-Phase Systems at
# T = 298.15 K. Journal of Chemical & Engineering Data, v. 55, n. 11, 
# p. 4741-4745, 2010/11/11 2010. ISSN 0021-9568. 
# (\href{https://pubs.acs.org/doi/full/10.1021/je100356s?src=recsys}
# {ACS Publications})
# 
# XIE, X.  et al. Liquidb-liquid equilibrium of aqueous two-phase systems of 
# PPG400 and biodegradable salts at temperatures of (298.15, 308.15, and 318.15)
#  K. Journal of Chemical & Engineering Data, v. 55, n. 8, p. 2857-2861,  2010. 
#  ISSN 0021-9568. 
# (\href{https://pubs.acs.org/doi/abs/10.1021/je901019t}{ACS Publications})
AQSys.default <- function(dataSET, modelName = "merchuk", Order="xy", ...) {
  # arrange data and guarantee R converted it to numbers but dont switch 
  # columns to prevent incompatibility with pre-existent functions 
  dataSET <- toNumeric(dataSET, Order)
  # each switch option calls a correspondent equation to fit dataSET
  # equations are functions declared in AQSysFormulas.R
  if (modelName %in% names(AQSysList(TRUE))) {
    ans <- do.call(modelName, list(dataSET))
  } else{
    AQSys.err("0")
  }
  # return fitting parameters ans statistical data
  return(ans)
}
#
###############################################################################
#' @rdname AQSys.data
#' @title Dataset and Fitted Function plot
#' @description The function returns a plot after fitting a dataset to the
#'  mathematical descriptor chosen by the user.
#' @details This version uses the plot function and return a regular orthogonal
#'  plot.
#' @method AQSys data
#' @export AQSys.data
#' @export
#' 
#' @param dataSET - Binodal Experimental data that will be used in the nonlinear
#'  fit. It might hold multiple systems stacked side-by-side. [type:data.frame]
#' @param modelName - Character String specifying the nonlinear empirical 
#' equation to fit data. [type:String]
#' The default method uses Merchuk's equation. Other mathematical descriptors 
#' can be listed using AQSysList().
#' @param Order Defines how the data is organized in the Worksheet. Use "xy" 
#' whether the first column corresponds to the lower phase fraction and "yx" 
#' whether the opposite. [type:String]
#' @param xmax Maximum value for the Horizontal axis' 
#' value - optional [type:double]
#' @param ymax Maximum value for the Vertical axis' 
#' value - optional [type:double]
#' @param ... Additional optional arguments. None are used at present.
#' @return return a data.frame with data fitted using the chosen mathematical 
#' descriptor.
#' @examples
#' # Populating variable dataSET with binodal data
#' dataSET <- peg4kslt[ , 1:2]
#' # Fitting dataSET using Merchuk's function
#' data <- AQSys.data(dataSET, Order = "xy")
AQSys.data <- function(dataSET,
                        modelName = "merchuk",
                        Order = "xy",
                        xmax = "",
                        ymax = "",
                        ...
) {
  # guarantee all lines are valid (non-na and numeric)
  dataSET <- toNumeric(dataSET, Order)
  # Calculate aesthetically better xmax and ymax
  if ((xmax == "") |
      ((xmax > 1) &
       (max(dataSET[1]) <= 1)) |
      ((xmax > 100) &
       (max(dataSET[1]) <= 100) &
       (max(dataSET[1]) >= 10)) |
      (xmax < round(max(dataSET[1]) / 0.92, 1))) {
    
    xmax <- ceiling(round(max(dataSET[1]) / 0.92, 1)/5)*5
  }
  if ((ymax == "") |
      ((ymax > 1) &
       (max(dataSET[2]) <= 1)) |
      ((ymax > 100) &
       (max(dataSET[2]) <= 100) &
       (max(dataSET[2]) >= 1)) |
      (ymax < round(max(dataSET[2]) / 0.92, 1))) {
    
    ymax <- ceiling(round(max(dataSET[2]) / 0.92, 1)/5)*5
  }
  # Select which model will be used to generate the plot. Function return list 
  # of plots and respective number of parameters
  models_npars <- AQSysList(TRUE)
  # choose model based in the user choice or standard value
  Fn <- ifelse(modelName %in% names(models_npars), 
               AQSys.mathDesc(modelName), AQSys.err("0"))
  # calculate coefficients using the non-linear regression
  CoefSET <- summary(AQSys(dataSET, modelName))$coefficients[, 1]
  #
  Xs <- as.data.frame(ifelse((xmax != max(dataSET[1])), 
                             rbind(dataSET[1], xmax), dataSET[1]))
  rawdt <- setNames(data.frame(Xs, Fn(CoefSET, Xs)), c("XC", "YC"))
  #
  return(rawdt)
}
###############################################################################
#' @rdname AQSys.plot
#' @title Dataset and Fitted Function plot
#' @description The function returns a plot after fitting a dataset to the 
#' mathematical descriptor chosen by the user.
#' @details This version uses the plot function and return a regular orthogonal 
#' plot.
#' @method AQSys plot
#' @export AQSys.plot
#' @export
#' @param dataSET - Binodal Experimental data that will be used in the nonlinear
#'  fit. It might hold multiple systems stacked side-by-side. [type:data.frame]
#' @param xlbl Plot's Horizontal axis label.
#' @param ylbl Plot's Vertical axis label.
#' @param main Legacy from plot package. For more details, 
#' see \code{\link{plot.default}}
#' @param col Legacy from plot package. For more details, 
#' see \code{\link{plot.default}}
#' @param type Legacy from plot package. For more details, 
#' see \code{\link{plot.default}}
#' @param cex Legacy from plot package. For more details, 
#' see \code{\link{plot.default}}
#' @param cexlab Legacy from plot package. For more details, 
#' see \code{\link{plot.default}}
#' @param cexaxis Legacy from plot package. For more details, 
#' see \code{\link{plot.default}}
#' @param cexmain Legacy from plot package. For more details, 
#' see \code{\link{plot.default}}
#' @param cexsub Legacy from plot package. For more details, 
#' see \code{\link{plot.default}}
#' @param modelName - Character String specifying the nonlinear empirical 
#' equation to fit data. The default method uses
#' Merchuk's equation. Other possibilities can be seen in AQSysList().
#' 
#' @param NP Number of points used to build the fitted curve. 
#' Default is 100. [type:Integer]
#' @param xmax Maximum value for the Horizontal axis' value
#' @param ymax Maximum value for the Vertical axis' value
#' @param Order Defines how the data is organized in the Worksheet. 
#' Use "xy" whether the first column corresponds to the lower phase fraction 
#' and "yx" whether the opposite.
#' @param save Save the generated plot in the disk using path and filename 
#' provided by the user. [type:Boulean]
#' @param HR Adjust Plot's text to be compatible with High Resolution 
#' size [type:Boulean]
#' @param filename Filename provided by the user to save a given 
#' plot. [type:String]
#' @param wdir The directory in which the plot file will be 
#' saved. [type:String]
#' @param silent save plot file without actually showing it to the 
#' user. [type:Boulean]
#' @param ... Additional optional arguments. None are used at present.
#' 
#' @return A plot containing the experimental data, the correspondent curve for
#'  the binodal in study and the curve's raw XY data.
#' @examples
#' #Populating variable dataSET with binodal data
#' dataSET <- peg4kslt[, 1:2]
#' # Plot dataSET using Merchuk's function
#' #
#' AQSys.plot(dataSET)
#' #
AQSys.plot <- function(dataSET,
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
                         modelName = "merchuk",
                         NP = 100,
                         xmax = "",
                         ymax = "",
                         Order = "xy",
                         save = FALSE,
                         HR = FALSE,
                         filename = NULL,
                         wdir = NULL,
                         silent = FALSE,
                         ...
)
  {
    #
    plot_image = NULL
    # verify and adjust data.frame names
    if (xlbl == "") {
      xlbl <- names(dataSET)[1]
    }
    if (ylbl == "") {
      ylbl <- names(dataSET)[2]
    }
    # guarantee all lines are valid (non-na and numeric)
    dataSET <- toNumeric(dataSET, Order)
    # Calculate aesthetically better xmax and ymax
    if ((xmax == "") |
        ((xmax > 1) &
         (max(dataSET[1]) <= 1)) |
        ((xmax > 100) &
         (max(dataSET[1]) <= 100) &
         (max(dataSET[1]) >= 10)) |
        (xmax < round(max(dataSET[1]) / 0.92, 1))) {
      
      xmax <- ceiling(round(max(dataSET[1]) / 0.92, 1)/5)*5
    }
    if ((ymax == "") |
        ((ymax > 1) &
         (max(dataSET[2]) <= 1)) |
        ((ymax > 100) &
         (max(dataSET[2]) <= 100) &
         (max(dataSET[2]) >= 1)) |
        (ymax < round(max(dataSET[2]) / 0.92, 1))) {
      
      ymax <- ceiling(round(max(dataSET[2]) / 0.92, 1)/5)*5
    }
    # If save=TRUE adjust variables and parameters to save plot
    wdir <- saveConfig(save, HR, filename, wdir, silent)
    # Select which model will be used to generate the plot. Function return 
    # list of plots and respective number of parameters
    models_npars <- AQSysList(TRUE)
    # choose model based in the user choice or standard value
    Fn <- ifelse(modelName %in% names(models_npars), AQSys.mathDesc(modelName),
                 AQSys.err("0"))
    # calculate coefficients using the non-linear regression
    CoefSET <- summary(AQSys(dataSET, modelName))$coefficients[, 1]
    # plot phase diagram using experimental data and with previously calculated
    #  parameters
    plot_image <- ggplot(data = dataSET, aes_string(x = "XC", y = "YC")) + 
      geom_point(shape = 8, size = 2) +
      theme_light() + xlab(paste(xlbl, "(%, m/m)")) + 
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
        limits = c(-0.5, ymax),
        breaks = seq(0, ymax, by = 5),
        labels = seq(0, ymax, by = 5)
      ) +
      scale_x_continuous(
        expand = c(0, 0),
        limits = c(0, xmax),
        breaks = seq(0, ymax, by = xmax / 10),
        labels = seq(0, ymax, by = xmax / 10)
      )
    #
    # add curve generated using regression parameters
    Xs <- as.data.frame(ifelse((xmax != max(dataSET[1])), 
                               rbind(dataSET[1], xmax), dataSET[1]))
    rawdt <- setNames(data.frame(Xs, Fn(CoefSET, Xs)), c("XC", "YC"))
    #
    plot_image <-
      plot_image + geom_line(
        data = rawdt,
        aes_string(x = "XC", y = "YC"),
        color = "red",
        linetype = "solid",
        size = 1
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
    # make available data from fitted curve to user. Function returns it
    #  silently but user can get data using simple assign '<-'
    if (silent == FALSE) {
      print(plot_image)
      invisible(rawdt)
    } else {
      invisible(plot_image)
    }
  }
###############################################################################
