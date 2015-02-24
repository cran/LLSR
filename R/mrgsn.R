#' @import rootSolve
require(rootSolve)
#'Murugesan & Perumalsamy (2005) Method
#' @rdname mrgsn
#' @name mrgsn
#' @description Murugesan and Perumalsamy published a paper in the Journal of Chemical Engineering Data describing
#'  how they used a set of equations to fit and study a PEG/Salt Aqueous Two-Phase System. They also correlated
#'   Tieline data using the Othmer-Tobias and Bancraft equations and the osmotic virial equation (the last one
#'    wasn't implemented yet). This R method implements the methodology described in the referred article.
#' @references Murugesan, T. and M. Perumalsamy, Liquid-Liquid Equilibria of Poly(ethylene glycol) 2000 + Sodium
#'  Citrate + Water at (25, 30, 35, 40, and 45) C. Journal of Chemical & Engineering Data, 2005. 50(4): p. 1392-1395.
#' @export
#' @seealso \itemize{
#'\item \code{\link{mrgsn.default}}
#'\item \code{\link{mrgsn.othmer}}
#'\item \code{\link{mrgsn.bancroft}}
#'\item \code{\link{mrgsn.plot}}
#' }
mrgsn <- function(...) UseMethod("mrgsn")
#'
#'
#' @rdname mrgsn
#' @title Marugesan's nonlinear Equation
#' @description Perform a nonlinear regression fit in order to determine the equation's parameters.
#' @details The function returns three parameters after fitting experimental data to the equation
#'  proposed Murugesan & Perumalsamy (2005)
#' @param XYdt - Binodal Experimental data that will be used in the nonlinear fit
#' @param ... Additional optional arguments. None are used at present.
#' @method mrgsn default
#' @export
#' @return Parameters P and Statistical data
#' @examples 
#' #Populating variable XYdt with binodal data
#' XYdt <- peg4kslt[,1:2] 
#' #Fitting XYdt using Murugesan's function
#' mrgsn(XYdt)
mrgsn.default <- function(XYdt,...){
  names(XYdt)<-c("XC","YC")
  controll<-nls.control(maxiter=50, tol=1e-10, minFactor = 1/1024, printEval = FALSE, warnOnly = FALSE)
  FFn <- nls(
    YC~ A+B*(XC)^0.5+C*XC,
    start=list(A=50,B=1,C=0),
    data=XYdt,na.exclude)    
  FFn
}
#' @rdname mrgsn.othmer
#' @title Othmer's Equation - Tieline's correlation
#' @description Othmer's equation to correlate tieline's data applying the lever's rule.
#' @references Othmer, D.F. and P.E. Tobias, Liquid -Liquid Extraction Data -Toluene and Acetaldehyde Systems.
#' Industrial & Engineering Chemistry, 1942. 34(6): p. 690-692.
#' @method mrgsn othmer
#' @param ... Additional optional arguments. None are used at present.
#' @param TLdt - Tieline Experimental data that will be used in the nonlinear fit
#' @export
#' @return Parameters K, n and Statistical data
#' @examples 
#' # TLdt is a data.frame which contains series of Tieline's mass fraction
#' # (upper-rich component, bottom-rich component and water)
#' # Each column in the data.frame represents a series of one component mass fraction
#' # For example, an empty data.frame for four tielines can be obtaining using:
#' TLdt<-matrix(NA,nrow=4,ncol=6)
#' # Variables order must follows the sequence presented below:
#' # "mfXt","mfYt","mfXb","mfYb","mfWt","mfWb"
#' # In which: mf stands for mass fraction; X and Y for the component
#' # rich in bottom and upper phase, respectively; t or b for top and
#' # bottom phases and W for water.
#' # Then you just need to load the data.frame in the function:
#' \dontrun{
#' mrgsn.othmer(TLdt)
#'}
mrgsn.othmer <- function(TLdt,...){
  #
  TLdt<-as.data.frame(TLdt)
  names(TLdt)<-c("mfXt","mfYt","mfXb","mfYb","mfWt","mfWb")
  #
  suppressWarnings(FFn <- nls(
    log((1-mfYt)/mfYt) ~ log(K*(((1-mfXb)/mfXb))^n),
    start=list(n=1,K=1), 
    algorithm="port",
    lower=10^-10,
    data=TLdt,na.exclude))
  FFn
}
#' @name mrgsn.bancroft
#' @title Bancroft's Potential Equation - Tieline's correlation
#' @description Bancroft's equation to correlate tieline's data.
#' @references Othmer, D.F. and P.E. Tobias, Liquid-Liquid Extraction Data -Toluene and Acetaldehyde Systems.
#'  Industrial & Engineering Chemistry, 1942. 34(6): p. 690-692.
#' @method mrgsn bancroft
#' @export
#' @param ... Additional optional arguments. None are used at present.
#' @param TLdt - Tieline Experimental data that will be used in the nonlinear fit
#' @return Parameters K1, r and Statistical data
#' @examples 
#' # TLdt is a data.frame which contains series of Tieline's mass fraction
#' # (upper-rich component, bottom-rich component and water)
#' # Each column in the data.frame represents a series of one component mass fraction
#' # For example, an empty data.frame for four tielines can be obtaining using:
#' TLdt<-matrix(NA,nrow=4,ncol=6)
#' # Variables order must follows the sequence presented below:
#' # "mfXt","mfYt","mfXb","mfYb","mfWt","mfWb"
#' # In which: mf stands for mass fraction; X and Y for the component
#' # rich in bottom and upper phase, respectively; t or b for top and
#' # bottom phases and W for water.
#' # Then you just need to load the data.frame in the function:
#' \dontrun{
#' mrgsn.bancroft(TLdt)
#'}
mrgsn.bancroft <- function(TLdt,...){
  #
  TLdt<-as.data.frame(TLdt)
  names(TLdt)<-c("mfXt","mfYt","mfXb","mfYb","mfWt","mfWb")
  #
  #suppressWarnings(
    FFn <- nls(
    log((mfWb/mfXb)) ~ log(K1*((mfWt/mfYt)^r)),
    start=list(r=1,K1=1),
    algorithm="port",
    lower=10^-2,
    data=TLdt,na.exclude)
  #)
  FFn
}
#
#' @name mrgsn.plot
#' @title This functions plot the Dataset and its Fitted curve
#' @description The function returns a plot after fitting a dataset to the equation described by Murugesan.
#' @details This version uses the plot function and return a regular bidimensional plot. Future versions will
#'  include a ternary diagram and more formal formatting.
#' @method mrgsn plot
#' @export
#' @param ... Additional optional arguments. None are used at present.
#' @param XYdt - Standard bidimensional data.frame used in most of functions available in this package. [type::data.frame]
#' @param xlbl - Plot's Horizontal axis label. If not set, It will admit the system under study is a PEG-Salt System.   [type:string]
#' @param ylbl - Plot's Vertical axis label. If not set, It will admit the system under study is a PEG-Salt System.   [type:string]
#' @param main - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param col - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param type - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cex - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexlab - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexaxis - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexmain - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param cexsub - Legacy from plot package. For more details, see \code{\link{plot.default}}
#' @param xmax - Maximum value for the Horizontal axis' value  [type:double]
#' @param ymax - Maximum value for the Vertical axis' value  [type:double]
#' @return A plot containing the experimental data and the correspondent curve for the binodal in study.
#' @examples 
#' #Populating variable XYdt with binodal data
#' XYdt <- peg4kslt[,1:2] 
#' #Plot XYdt using Murugesan's function
#' \dontrun{
#' mrgsn.plot(XYdt)
#' }
mrgsn.plot <- function  (XYdt, xlbl="Salt Fraction (w/w)", ylbl="PEG Fraction (w/w)", main="Title", col="blue", type="p",
                         cex=2.5, cexlab=2.5, cexaxis=2.5, cexmain=2.5, cexsub=2.5,xmax=0.4,ymax=0.5,...)
  {
  #
  par(mar = c(6,6,6,4) + 0.1)
  #
  plot(XYdt, xlab=xlbl,ylab=ylbl,main=main,col=col, type=type,
       cex=cex,cex.lab=cexlab, cex.axis=cexaxis, cex.main=cexmain,
       cex.sub=cexsub,xlim=c(0,xmax),ylim=c(0,ymax))  
  #
  Smmry<-summary(mrgsn(XYdt))
  #
  x<-sort(runif(500,0.001,xmax))
  #
  Fn <- function(A,B,C,XC){
    A+B*(XC)^0.5+C*XC
  }
  #
  curve(Fn(Smmry$coefficients[1],
                Smmry$coefficients[2],
                Smmry$coefficients[3],x),
        add=TRUE)
}