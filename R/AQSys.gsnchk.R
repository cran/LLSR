#' @import rootSolve
#' @rdname AQSys.gsnchk
#' @name AQSys.gsnchk
#' @title Merchuk's Method of calculating phase's composition in a tie-line applied to Murugesan's bimodal function
#' @description Merchuk et al. described a very straightforward method to calculate the concentration of each
#'component in the tie-line giving only its global composition and phase's properties (such as volume and density).
#'However, other researchers relate to have achieved better fitting results using Murugesan's functions.
#'This method offers Merchuk's methodology to calculate tie-line's composition but using the equation proposed by Murugesan.
# Merchuk's method of calculating phase's composition in a tie-line applied to Murugesan's bimodal function
#' @details Using the bimodal data, the global composition of a chosen tieline and its phases properties (more precisely
#'each phase density and volume). Using the data included in LLSR package the function couldn't achieve steady-state and
#'consecutively have a poor convergence tolerance. Use for your own risk.
#' @references 
#' MURUGESAN, T.; PERUMALSAMY, M. Liquid-Liquid Equilibria of Poly(ethylene glycol) 2000 + Sodium Citrate + Water at (25, 30, 35, 40, and 45) C. Journal of Chemical & Engineering Data, v. 50, n. 4, p. 1392-1395, 2005/07/01 2005. ISSN 0021-9568. 
#' (\href{https://www.doi.org/10.1021/je050081k}{ACS Publications})
#' 
#' MERCHUK, J. C.; ANDREWS, B. A.; ASENJO, J. A. Aqueous two-phase systems for protein separation: Studies on phase inversion. Journal of Chromatography B: Biomedical Sciences and Applications, v. 711, n. 1-2, p. 285-293,  1998. ISSN 0378-4347.
#' (\href{https://www.doi.org/10.1016/s0378-4347(97)00594-x}{ScienceDIrect})
#' @export AQSys.gsnchk
#' @export
#' @param dataSET - Standard bidimensional data.frame used in most of functions available in this package.
#' @param Xm - Component X's concentration in the tieline's global composition.
#' @param Ym - Component Y's concentration in the tieline's global composition.
#' @param Vt - Tieline's TOP phase volume.
#' @param Vb - Tieline's BOTTOM phase volume.
#' @param dyt - Tieline's TOP phase density
#' @param dyb - Tieline's BOTTOM phase density
#' @param ... Additional optional arguments. None are used at present.
#' @return sysres - The function returns the Critical Point (X,Y), Tieline Length (TLL), Tieline's Equivolume point (xVRe2o,yVRe2o),
#'and Tieline's Slope.
#' @examples
#' #
#' dataSET <- peg4kslt[,1:2]
#' #
#' Xm <- peg4kslt[2,3]
#' Ym <- peg4kslt[2,4]
#' Vt <- peg4kslt[2,5]
#' Vb <- peg4kslt[2,6]
#' dyt <- peg4kslt[2,7]
#' dyb <- peg4kslt[2,8]
#' #
#' AQSys.gsnchk(dataSET,Xm,Ym,Vt,Vb,dyt,dyb)
AQSys.gsnchk <- function(dataSET,Xm,Ym,Vt,Vb,dyt,dyb,...) {
  # Fit dataSET data to Murugesan's equation and store it in Smmry
  Smmry <- summary(murugesan(dataSET))
  # extract regression parameters from Smmry
  A <- Smmry$coefficients[1]
  B <- Smmry$coefficients[2]
  C <- Smmry$coefficients[3]
  # Calculate alfa for a given system composition
  alfa <- Vt * dyt / (Vt * dyt + Vb * dyb)
  # the system of equations below was uses the method described by Merchuk
  # in the manuscript Merchuk, J.C., B.A. Andrews, and J.A. Asenjo. (1998),
  # Aqueous two-phase systems for protein separation: Studies on phase inversion.
  # Journal of Chromatography B: Biomedical Sciences and Applications. 711(1-2):
  #p. 285-293. to calculate tieline composition using murugesan's equation.
  # the lines below set the equation's system
  sys <- function(x) {
    F1 <- A + B * (x[2]) ^ 0.5 + C * x[2] - x[1]
    F2 <- A + B * (x[4]) ^ 0.5 + C * x[4] - x[3]
    F3 <- (Ym / alfa) - ((1 - alfa) / alfa) * x[3] - x[1]
    F4 <- (Xm / alfa) - ((1 - alfa) / alfa) * x[4] - x[2]
    #
    c(F1 = F1, F2 = F2, F3 = F3, F4 = F4)
  }
  # solve the system of equation for a given set of guess and restricting of positive
  # only results
  (sysres <- multiroot(f = sys, start = c(1,0,0,1), positive = TRUE))
  # Calculate the tieline length and store it in sysres under the TLL alias
  sysres$TLL <- sqrt((sysres$root[1] - sysres$root[3]) ^ 2 + (sysres$root[2] - sysres$root[4]) ^ 2)
  # set alfa to 0.5 to calculate concentration at equivolume point
  alfaVRe2o <- 0.5
  # calculate the system composition at equivolume
  sysres$yVRe2o <- alfaVRe2o * (sysres$root[1] + sysres$root[3] * ((1 - alfaVRe2o) / alfaVRe2o))
  sysres$xVRe2o <- alfaVRe2o * (sysres$root[2] + sysres$root[4] * ((1 - alfaVRe2o) / alfaVRe2o))
  # set var name for root results (phase's composition for a given tieline)
  names(sysres$root) <- c("YT","XT","YB","XB")
  # calculate and store tieline's slope
  sysres$S <- (sysres$root["YT"] - sysres$root["YB"]) / (sysres$root["XT"] - sysres$root["XB"])
  # removing Slope's header to make easier its retrieve
  names(sysres$S) <- NULL
  # return all calculated parameters
  sysres
}
