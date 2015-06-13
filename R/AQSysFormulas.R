require(rootSolve)
#
mrchk <- function(XYdt,P1=10,P2=1,P3=0,...){
  #
  names(XYdt)<-c("XC","YC")
  #
  FFn <- nls(
    YC~ P1*exp(P2*(XC^(0.5))-P3*(XC^3)),
    start=list(P1=P1,P2=P2,P3=P3),
    data=XYdt,na.exclude)    
  FFn
}

mrgsn <- function(XYdt,...){
  #
  names(XYdt)<-c("XC","YC")
  ##  controll<-nls.control(maxiter=50,
  ##  tol=1e-10, minFactor = 1/1024,
  ##  printEval = FALSE, warnOnly = FALSE)
  FFn <- nls(
    YC~ P1+P2*(XC)^0.5+P3*XC,
    start=list(P1=50,P2=1,P3=0),
    data=XYdt,na.exclude)    
  FFn
}

tello <- function(XYdt,...){
  #
  df.sys <- diff(XYdt[,1])/diff(XYdt[,2])
  nrow.sys <- nrow(XYdt)
  S<-smooth.spline(df.sys)$fit$coef[1:nrow.sys]
  XC <- XYdt[,1]
  #
  coef.est <- LLSRxy(XC,S)
  names(coef.est)<-c("XC","S")
  #
  FFnEst <- nls(
    S ~ (P2/P1)+XC/P1,
    start=list(P1=-.1,P2=.001),
    data=coef.est,na.exclude)
  #
  coefEst <- summary(FFnEst)
  coef.1 <- coefEst$coefficients[1]
  coef.2 <- coefEst$coefficients[2]
  #
  x <- mean(XYdt[,1])
  y <- mean(XYdt[,2])
  coef.3 <- log(exp(y)*((x+coef.2)^(-coef.1)))
  #
  #if (coef.2 < 0) coef.2 <- -coef.2
  if (coef.2 < -min(XYdt[,1])) coef.2 <- -(min(XYdt[,1])-0.001)
  #
  #cat(coef.1," ",coef.2," ",coef.3," ")
  #
  names(XYdt)<-c("XC","YC")
  FFn <- nls(
    #XC ~ exp((YC - P3)/P1) - P2,
    YC ~ P1*log(XC + P2) + P3,
    start = list(P1 = coef.1, P2 = coef.2, P3 = coef.3),
    #algorithm = "port",
    #lower = c(P1 = -Inf, P2 = -(min(XYdt[,1])-0.001), P3 = -Inf),
    data = XYdt, na.action = na.exclude) 
  FFn
}
