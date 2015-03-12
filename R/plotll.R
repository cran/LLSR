#' @rdname plotll
#' @title This functions plot a curve based in the chosen model and its parameters.
#' @description The function returns a plot after using the parameters and model given by the user.
#' @details The function owns predefined set of equations that can be seen below and must be used, with adequated parameters,
#' to return a plot which represent the chosen model.
#' @export plotll
#' @param model - Model chosen to be used: merchuk, murugesan [type:string]
#' @param param - Model's parameters [type::data.frame]
#' @param xlbl - Plot's Horizontal axis label. If not set, It will admit the system under study is a PEG-Salt System.  [type:string]
#' @param ylbl - Plot's Vertical axis label. If not set, It will admit the system under study is a PEG-Salt System.  [type:string]
#' @param main - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param col - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param type - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cex - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexlab - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexaxis - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexmain - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexsub - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param xmax - Maximum value for the Horizontal axis' value (bottom-rich component)  [type:double]
#' @param mpl - You can overlay curves in a single plot. Default is FALSE. [type::LOGIC]
#' @param HR - Magnify Plot's text to be compatible with High Resolution size [type:Boulean]
#' @param NP - Number of points used to build the fitted curve. Default is 100. [type:Integer]
#' @return A plot using the input model within the chosen interval and the curve's raw XY data. If no interval is selected, xmax=0.4.
#' @examples
#' \dontrun{
#' plotll("murugesan",as.data.frame(c(0.90,-3.48,2.92)),mpl=TRUE,col="red")
#' }
plotll <- function  (model,param, xlbl = "Salt Fraction (w/w)", ylbl = "PEG Fraction (w/w)", main = "Title", col = "black", type = "p",
                     cex = 1, cexlab = 1, cexaxis = 1, cexmain = 1, cexsub = 1, xmax = 0.4, mpl = FALSE, HR = FALSE, NP = 100)
{
  #
  if (HR==TRUE){
    par(mar = c(6,6,6,4) + 0.1)
    cex=2.5
    cexlab=2.5
    cexaxis=2.5
    cexmain=2.5
    cexsub=2.5
  }else{
    par(mar = c(5, 4, 4, 2) + 0.1)
    cex=1
    cexlab=1
    cexaxis=1
    cexmain=1
    cexsub=1
  }
  #
  x<-sort(runif(NP,0.001,xmax))
  if (model=="merchuk"){
    Fn1 <- function(P1,P2,P3,x){
      P1*exp(P2*(x^(0.5))-P3*(x^3))
    }
    rawdt <- curve(Fn1(param[1,1], param[2,1], param[3,1],x),col=col,add=mpl,xlim=c(0,xmax),xlab=xlbl,ylab=ylbl)
  } else{
    Fn2 <- function(A,B,C,XC){
      A+B*(XC)^0.5+C*XC
    }
    rawdt <- curve(Fn2(param[1,1], param[2,1], param[3,1],x),col=col,add=mpl,xlim=c(0,xmax),xlab=xlbl,ylab=ylbl)
  }
  names(rawdt)<-c("XC","YC")
  rawdt<-as.data.frame(rawdt)
  invisible(rawdt)
}