#'Merchuk's Equation to fit Binodal Experimental Data
#' @rdname mrchk
#' @name mrchk
#' @description Merchuk et al. published a paper in the Journal of 'Chromatography B: Biomedical Sciences and Applications' 
#' relating the use of a nonlinear equation to fit and study phase inversion and protein separation using a PEG400/Phosphate 
#' Aqueous Two-Phase System. They also calculated Tieline composition by using a set of four equations and the relationship 
#' of Tieline's Slope with tielines parameters. This R method implements the methodology described in the referred article.
#' @references Merchuk, J.C., B.A. Andrews, and J.A. Asenjo, Aqueous two-phase systems for protein separation: Studies on 
#' phase inversion. Journal of Chromatography B: Biomedical Sciences and Applications, 1998. 711(1-2): p. 285-293.
#' @return Parameters P and Statistical data
#' @export
#' @seealso \itemize{
#' \item \code{\link{mrchk.default}}
#' \item \code{\link{mrchk.plot}}
#' \item \code{\link{mrchk.tielines}} 
#' \item \code{\link{mrchk.crpt}}
#' \item \code{\link{gsnchk}}
#' }
mrchk <- function(XYdt,...) UseMethod("mrchk")
#' @rdname mrchk
#' @title Merchuk's nonlinear Equation
#' @description Perform a nonlinear regression fit in order to determine the equation's parameters.
#' @details The function returns three parameters after fitting experimental data to the equation proposed by Merchuk et al.
#' @param XYdt - Binodal Experimental data that will be used in the nonlinear fit
#' @param ... Additional optional arguments. None are used at present.
#' @param P1,P2,P3 Merchuk's Parameters
#' @method mrchk default
#' @export
#' @return Parameters P and Statistical data
#' @examples
#' #Populating variable XYdt with binodal data
#' XYdt <- peg4kslt[,1:2] 
#' #Fitting XYdt using Murugesan's function
#' mrchk(XYdt)
mrchk.default <- function(XYdt,P1=10,P2=1,P3=0,...){
  #
  names(XYdt)<-c("XC","YC")
  #
  FFn <- nls(
    YC~ P1*exp(P2*(XC^(0.5))-P3*(XC^3)),
    start=list(P1=P1,P2=P2,P3=P3),
    data=XYdt,na.exclude)    
  FFn
}
# MERCHUK PLOT TEST FUNCTION
#' @rdname mrchk.plot
#' @title Dataset and Fitted Function plot
#' @description The function returns a plot after fitting a dataset to the equation described by Merchuk.
#' @details This version uses the plot function and return a regular bidimensional plot. Future versions will include a ternary 
#' diagram and more formal formatting.
#' @method mrchk plot
#' @export
#' @param ... Additional optional arguments. None are used at present.
#' @param XYdt - Binodal Experimental data that will be used in the nonlinear fit
#' @param xlbl = Plot's Horizontal axis label. If not set, It will admit the system under study is a PEG-Salt System.
#' @param ylbl = Plot's Vertical axis label. If not set, It will admit the system under study is a PEG-Salt System.
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
#' @return A plot containing the experimental data and the correspondent curve for the binodal in study.
#' @examples
#' #Populating variable XYdt with binodal data
#' XYdt <- peg4kslt[,1:2] 
#' #Plot XYdt using Merchuk's function
#' \dontrun{
#' mrchk.plot(XYdt)
#' }
mrchk.plot <- function  (XYdt, xlbl="Salt Fraction (w/w)", ylbl="PEG Fraction (w/w)", main=NULL, col="blue", type="p",
                         cex=2.5, cexlab=2.5, cexaxis=2.5, cexmain=2.5, cexsub=2.5,xmax=0.4,ymax=0.5,...)
{
    #
    par(mar = c(6,6,6,4) + 0.1)
    #
    plot(XYdt, xlab=xlbl,ylab=ylbl,main=main,col=col, type=type,
         cex=cex,cex.lab=cexlab, cex.axis=cexaxis, cex.main=cexmain,
         cex.sub=cexsub,xlim=c(0,xmax),ylim=c(0,ymax))  
    #Test if summary exists
    Smmry<-summary(mrchk(XYdt))
    #
    x<-sort(runif(500,0.001,xmax))
    #
    Fn <- function(P1,P2,P3,SALT){
      P1*exp(P2*(SALT^(0.5))-P3*(SALT^3))
    }
    #
    curve(Fn(Smmry$coefficients[1],
                    Smmry$coefficients[2],
                    Smmry$coefficients[3],x),
          add=TRUE)
}
#' @rdname mrchk.tielines
#' @title Merchuk's Method - Tieline's Composition Calculation
#' @description Merchuk et al. described a very straightforward method to calculate the concentration of each component in the 
#' tieline giving only its global composition and phase's properties (such as volume and density). 
#' @details Using the binodal data, the global composition of a chosen tieline and its phases properties (more precisely each 
#' phase density and volume)
#' @method mrchk tielines
#' @export
#' @param XYdt - Binodal Experimental data that will be used in the nonlinear fit
#' @param Xm - Component X's concentration in the tieline's global composition.
#' @param Ym - Component Y's concentration in the tieline's global composition.
#' @param Vt - Tieline's TOP phase volume.
#' @param Vb - Tieline's BOTTOM phase volume.
#' @param dyt - Tieline's TOP phase density
#' @param dyb - Tieline's BOTTOM phase density
#' @param ... Additional optional arguments. None are used at present.
#' @return The function returns the Critical Point (X,Y), Tieline Length (TLL), Tieline's Equivolume point (xVRe2o,yVRe2o), 
#' and Tieline's Slope.
#' @examples
#' \dontrun{
#' mrchk.tielines(XYdt,Xm,Ym,Vt,Vb,dyt,dyb)
#' }
mrchk.tielines<- function(XYdt,Xm,Ym,Vt,Vb,dyt,dyb,...){
  #
  Smmry<-summary(mrchk(XYdt))
  #
  P1<-Smmry$coefficients[1]
  P2<-Smmry$coefficients[2]
  P3<-Smmry$coefficients[3]
  #
  alfa<-Vt*dyt/(Vt*dyt+Vb*dyb)
  #
  sys <- function(x) {
    F1 <- P1*exp(P2*x[2]^0.5-P3*x[2]^3) - x[1] 
    F2 <- P1*exp(P2*x[4]^0.5-P3*x[4]^3) - x[3]
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
#' @rdname mrchk.crpt
#' @title Merchuk's Method - Critical Point Calculation
#' @description Merchuk et al. described a very straightforward method to calculate the critical composition of a given Binodal 
#' curve and its Tielines.
#' @details Using the binodal data, tieline's Slopes (S), the composition of bottom-rich component in the bottom phase (XB) 
#' and an equation which stablish a relatioship between them, this function returns the critical composition of the binodal 
#' under study. When used within a iterative function, mrchk.tielines() can be used to obtain TLL and S and therefore 
#' calculate the critical composition.
#' @method mrchk crpt
#' @export
#' @param tldata - A data.frame with two columns containing a set of Tieline's Slopes (S)
#' and its bottom-rich component composition in the bottom phase (XB).
#' @param XYdt - Binodal Experimental data that will be used in the nonlinear fit
#' @param ... Additional optional arguments. None are used at present.
#' @return (XCrit,YCrit) - The function returns Tieline's Critical Point Composition
#' @examples
#' \dontrun{
#' mrchk.crpt(XYdt,tldata)
#' } 
mrchk.crpt <- function(XYdt,tldata,...){
  Smmry<-summary(mrchk(XYdt))
  #
  tldata<-as.data.frame(tldata)
  resTLLfit<-lm(S ~ poly(XB, 2, raw=TRUE),data=tldata)
  #
  fitC <- function(x) {
    F1= -x[1] + Smmry$coefficients[1]*exp(Smmry$coefficients[2]*x[2]^0.5-Smmry$coefficients[3]*x[2]^3)*((Smmry$coefficients[2]/(2*x[2]^0.5))-3*Smmry$coefficients[3]*x[2]^2)
    F2= -x[1] + resTLLfit$coefficients[1]+resTLLfit$coefficients[2]*x[2]+resTLLfit$coefficients[3]*x[2]^2
    c(F1=F1,F2=F2)
  }
  # 
  (sres <- multiroot(f = fitC, start = c(.1,.1)))
  XCrit<-sres$root[2]
  YCrit<-Smmry$coefficients[1]*exp(Smmry$coefficients[2]*(XCrit^(0.5))-Smmry$coefficients[3]*(XCrit^3))
  #
  sres$XCrit<-XCrit
  sres$YCrit<-YCrit
  #
  sres
}