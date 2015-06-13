#' @rdname AQSys.tielines
#' @title Merchuk's Method - Tieline's Composition Calculation
#' @description Merchuk et al. described a very straightforward method to calculate the concentration of each component in the 
#' tieline giving only its global composition and phase's properties (such as volume and density). 
#' @details Using the binodal data, the global composition of a chosen tieline and its phases properties (more precisely each 
#' phase density and volume)
#' @method AQSys tielines
#' @export AQSys.tielines
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
#' AQSys.tielines(XYdt,Xm,Ym,Vt,Vb,dyt,dyb)
#' }
AQSys.tielines<- function(XYdt,Xm,Ym,Vt,Vb,dyt,dyb,...){
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
#' @rdname AQSys.crpt
#' @title Merchuk's Method - Critical Point Calculation
#' @description Merchuk et al. described a very straightforward method to calculate the critical composition of a given Binodal 
#' curve and its Tielines.
#' @details Using the binodal data, tieline's Slopes (S), the composition of bottom-rich component in the bottom phase (XB) 
#' and an equation which stablish a relatioship between them, this function returns the critical composition of the binodal 
#' under study. When used within a iterative function, mrchk.tielines() can be used to obtain TLL and S and therefore 
#' calculate the critical composition.
#' @method AQSys crpt
#' @export AQSys.crpt
#' @export
#' @param tldata - A data.frame with two columns containing a set of Tieline's Slopes (S)
#' and its bottom-rich component composition in the bottom phase (XB).
#' @param XYdt - Binodal Experimental data that will be used in the nonlinear fit
#' @param ... Additional optional arguments. None are used at present.
#' @return (XCrit,YCrit) - The function returns Tieline's Critical Point Composition
#' @examples
#' \dontrun{
#' AQSys.crpt(XYdt,tldata)
#' } 
AQSys.crpt <- function(XYdt,tldata,...){
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