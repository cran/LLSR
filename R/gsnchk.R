#' @import rootSolve
require(rootSolve)
#'Merchuk's Method of calculating phase's composition in a Tieline applied to Murugesan Binodal Function
#' @rdname gsnchk
#' @name gsnchk
#' @title Merchuk's Method of calculating phase's composition in a Tieline applied to Murugesan Binodal Function
#' @description Merchuk et al. described a very straightforward method to calculate the concentration of each
#'component in the tieline giving only its global composition and phase's properties (such as volume and density).
#'However, other researchers relate to have achieved better fitting results using Murugesan's functions.
#'This method offers Merchul's ways o calculate Tieline's composition but using the equation proposed by Murugesan. 
#' @details Using the binodal data, the global composition of a chosen tieline and its phases properties (more precisely
#'each phase density and volume). Using the data included in LLSR package the function couldn't achieve steady-state and
#'consecutively have a poor convergence tolerance. Use for your own risk.
#' @export gsnchk
#' @param XYdt - Standard bidimensional data.frame used in most of functions available in this package.
#' @param Xm - Component X's concentration in the tieline's global composition.
#' @param Ym - Component Y's concentration in the tieline's global composition.
#' @param Vt - Tieline's TOP phase volume.
#' @param Vb - Tieline's BOTTOM phase volume.
#' @param dyt - Tieline's TOP phase density
#' @param dyb - Tieline's BOTTOM phase density
#' @return sysres - The function returns the Critical Point (X,Y), Tieline Length (TLL), Tieline's Equivolume point (xVRe2o,yVRe2o),
#'and Tieline's Slope.
#' @examples
#' #
#' XYdt <- peg4kslt[,1:2]
#' #
#' Xm <- peg4kslt[2,3]
#' Ym <- peg4kslt[2,4]
#' Vt <- peg4kslt[2,5]
#' Vb <- peg4kslt[2,6]
#' dyt <- peg4kslt[2,7]
#' dyb <- peg4kslt[2,8]
#' #
#' gsnchk(XYdt,Xm,Ym,Vt,Vb,dyt,dyb)
gsnchk<- function(XYdt,Xm,Ym,Vt,Vb,dyt,dyb){
  #
  Smmry<-summary(mrgsn(XYdt))
  #
  A<-Smmry$coefficients[1]
  B<-Smmry$coefficients[2]
  C<-Smmry$coefficients[3]
  #
  alfa<-Vt*dyt/(Vt*dyt+Vb*dyb)
  #
  sys <- function(x) {
    F1 <- A+B*(x[2])^0.5+C*x[2] - x[1]
    F2 <- A+B*(x[4])^0.5+C*x[4] - x[3]
    F3 <- (Ym/alfa)-((1-alfa)/alfa)*x[3] - x[1]
    F4 <- (Xm/alfa)-((1-alfa)/alfa)*x[4] - x[2]
    #
    c(F1 = F1, F2 = F2, F3 = F3, F4 = F4)
  }
  #
  (sysres <- multiroot(f = sys, start = c(1,0,0,1),positive=TRUE))
  #
  sysres$TLL<-sqrt((sysres$root[1]-sysres$root[3])^2+(sysres$root[2]-sysres$root[4])^2)
  #
  alfaVRe2o <- 0.5
  #
  sysres$yVRe2o <- alfaVRe2o*(sysres$root[1]+sysres$root[3]*((1-alfaVRe2o)/alfaVRe2o))
  sysres$xVRe2o <- alfaVRe2o*(sysres$root[2]+sysres$root[4]*((1-alfaVRe2o)/alfaVRe2o))
  #
  names(sysres$root)<-c("YT","XT","YB","XB")
  sysres$S <- (sysres$root["YT"]-sysres$root["YB"])/(sysres$root["XT"]-sysres$root["XB"])
  names(sysres$S)<-NULL
  sysres
}

#' Dataset of experimental binodal data of an ATPS
#'
#' A dataset containing the experimental binodal data for a PEG/SALT Aqueous Two-Phases System (ATPS)
#' @name peg4kslt
#' @format A data.frame with 116 rows and 2 variables:
#' \describe{
#'   \item{XC}{Fraction of Salt}
#'   \item{YC}{Fraction of PEG}
#'   ...
#' }
#' @source \url{http://diegofcoelho.github.io/code/}
#' @examples
#' #
#' XYdt <- peg4kslt[,1:2]
#' #
#' Xdt<-peg4kslt[,1]
#' #
#' Ydt<-peg4kslt[,1]
NULL
#'
#' @rdname xyJn
#' @title Dataset aggregation Function
#' @description This function returns a data.frame from given two unidimensionals vectors, XC and YX.
#' @details The function returns a data.frame after merging two unidimensional vectors.
#' @export xyJn
#' @param ... Additional parameters - Not implemented
#' @param YC - Component Y's concentration in the TOP Y-rich phase.
#' @param XC - Component X's concentration in the BOTTOM X-rich phase.
#' @return XYdt - Standard bidimensional data.frame used in most of functions available in this package.
#' @examples
#' #
#' Xdt<-peg4kslt[,1]
#' #
#' Ydt<-peg4kslt[,1]
#' #
#' xyJn(Xdt,Ydt)
xyJn <- function(XC,YC,...){
  xc <- as.vector(XC)
  yc <- as.vector(YC)
  XYdt<-data.frame(XC=xc,YC=yc)
  XYdt
}
#' @rdname plotll
#' @title This functions plot a curve based in the chosen model and its parameters.
#' @description The function returns a plot after using the parameters ans model given by the user.
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
#' @return A plot using the input model within the chosen interval. If no interval is selected, xmax=0.4.
#' @examples
#' \dontrun{
#' plotll("merugesan",as.data.frame(c(0.90,-3.48,2.92)),mpl=TRUE,col="red")
#' }
plotll <- function  (model,param, xlbl="Salt Fraction (w/w)", ylbl="PEG Fraction (w/w)", main="Title", col="black", type="p",
                     cex=2.5, cexlab=2.5, cexaxis=2.5, cexmain=2.5, cexsub=2.5,xmax=0.4,mpl=FALSE)
{
  x<-sort(runif(500,0.001,xmax))
  if (model=="merchuk"){
      Fn1 <- function(P1,P2,P3,x){
        P1*exp(P2*(x^(0.5))-P3*(x^3))
      }
      curve(Fn1(param[1,1], param[2,1], param[3,1],x),col=col,add=mpl,xlim=c(0,xmax))
  } else{
      Fn2 <- function(A,B,C,XC){
        A+B*(XC)^0.5+C*XC
      }
      curve(Fn2(param[1,1], param[2,1], param[3,1],x),col=col,add=mpl,xlim=c(0,xmax))
  }
}