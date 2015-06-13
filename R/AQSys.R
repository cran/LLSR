####################################################################################################################
#' @import rootSolve
####################################################################################################################
require(rootSolve)
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
AQSys <- function(XYdt,...) UseMethod("AQSys")
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
AQSys.default <- function(XYdt, mathDesc="merchuk",...){
  switch(mathDesc,
         merchuk={
           ans<-mrchk(XYdt)
         },
         murugesan={
           ans<-mrgsn(XYdt)
         },
         tello={
           ans<-tello(XYdt)
         },
         AQSys.err("0")
  )
  return(ans)
}
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
#' @param XYdt - Binodal Experimental data that will be used in the nonlinear fit
#' @param xlbl = Plot's Horizontal axis label.
#' @param ylbl = Plot's Vertical axis label. 
#' @param main - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param col - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param type - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cex - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexlab - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexaxis - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexmain - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexsub - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param xmax - Maximum value for the Horizontal axis' value
#' @param ymax - Maximum value for the Vertical axis' value
#' @param HR - Magnify Plot's text to be compatible with High Resolution size [type:Boulean]
#' @param NP - Number of points used to build the fitted curve. Default is 100. [type:Integer]
#' @return A plot containing the experimental data, the correspondent curve for the binodal in study and the curve's raw XY data.
#' @examples
#' #Populating variable XYdt with binodal data
#' XYdt <- peg4kslt[,1:2] 
#' #Plot XYdt using Merchuk's function
#' #
#' AQSys.plot(XYdt)
#' #
AQSys.plot <- function  (XYdt, xlbl = "", ylbl = "", main = NULL, col = "blue", type = "p",
                         cex = 1, cexlab = 1, cexaxis = 1, cexmain = 1, cexsub = 1,
                         xmax = 0.4, ymax = 0.5, HR = FALSE, NP = 100, mathDesc = "merchuk",...)
{
  #
  AQSysHR(TRUE)
  #
  switch(mathDesc,
         merchuk={
           CoefSET <- summary(mrchk(XYdt))$coefficients[, 1]
           Fn <- AQSys.mathDesc("merchuk")
         },
         murugesan={
           CoefSET <- summary(mrgsn(XYdt))$coefficients[, 1]
           Fn <- AQSys.mathDesc("murugesan")
         },
         tello={
           CoefSET <- summary(tello(XYdt))$coefficients[, 1]
           Fn <- AQSys.mathDesc("tello")
         },
        AQSys.err("0")
  )
  #
  plot(XYdt, xlab = xlbl, ylab = ylbl, main = main, col = col, type = type,
       cex = cex, cex.lab = cexlab, cex.axis = cexaxis, cex.main = cexmain,
       cex.sub = cexsub, xlim = c(0,xmax), ylim = c(0,ymax))  
  #
  x <- sort(runif(NP,0.001,xmax))
  #
  #SWITCH WORKS ONLY FOR THREE PARAMETER'S EQUATIONS.
  #IF NECESSARY A HIGHER NUMBER, INSERT CONDITIONAL BELOW.
  #Maybe change it to have as input the whole coefficient set?
  #a<-summary(mrchk(peg4kslt[,1:2]))$coefficients[,1]
  #
  rawdt <- curve(Fn(CoefSET,x),
                 add=TRUE, n = NP)
  names(rawdt)<-c("XC","YC")
  rawdt<-as.data.frame(rawdt)
  invisible(rawdt)
}
####################################################################################################################
