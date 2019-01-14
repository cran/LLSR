#' @rdname AQSys.tielines
#' @title Merchuk's Method - tie-line's Composition Calculation
#' @description Merchuk et al. described a very straightforward method to calculate the concentration of each component in the
#' tieline giving only its global composition and phase's properties (such as volume and density).
#' @details Using the binodal data, the global composition of a chosen tieline and its phases properties (more precisely each
#' phase density and volume)
#' @method AQSys tielines
#' @export AQSys.tielines
#' @export
#' @param dataSET - Binodal Experimental data that will be used in the nonlinear fit
#' @param Xm - Component X's concentration in the tieline's global composition.
#' @param Ym - Component Y's concentration in the tieline's global composition.
#' @param Vt - Tieline's TOP phase volume.
#' @param Vb - Tieline's BOTTOM phase volume.
#' @param dyt - Tieline's TOP phase density
#' @param dyb - Tieline's BOTTOM phase density
#' @param Wt - ATPS upper phase weight
#' @param Wb - ATPS bottom phase weight
#' @param byW - Use weight (TRUE) or volume and density (FALSE) during lever arm rule calculation.
#' @param ... Additional optional arguments. None are used at present.
#' @return The function returns the Critical Point (X,Y), Tieline Length (TLL), Tieline's Equivolume point (xVRe2o,yVRe2o),
#' and Tieline's Slope.
#' @examples
#' \dontrun{
#' AQSys.tielines(dataSET, Xm, Ym, Vt, Vb, dyt, dyb, wt, wb, byW = FALSE)
#' }
#' @references 
#' MERCHUK, J. C.; ANDREWS, B. A.; ASENJO, J. A. Aqueous two-phase systems for protein separation: Studies on phase inversion. Journal of Chromatography B: Biomedical Sciences and Applications, v. 711, n. 1-2, p. 285-293,  1998. ISSN 0378-4347.
#' (\href{https://www.doi.org/10.1016/s0378-4347(97)00594-x}{ScienceDIrect})
#' 
AQSys.tielines <-
  function(dataSET,
           Xm,
           Ym,
           Vt = NULL,
           Vb = NULL,
           dyt = NULL,
           dyb = NULL,
           Wt = NULL,
           Wb = NULL,
           byW = TRUE,
           ...) {
  # Fit dataSET data to Merchuk's equation and store it in Smmry
  Smmry <- summary(merchuk(dataSET))
  # extract regression parameters from Smmry
  P1 <- Smmry$coefficients[1]
  P2 <- Smmry$coefficients[2]
  P3 <- Smmry$coefficients[3]
  #
  if (byW & !is.null(c(Wt, Wb))) {
    # Calculate alfa for a given system composition
    alfa <- wt / (wt + wb)
  } else{
    AQSys.err("13")
  }
  #
  if (!byW & !is.null(c(Vt, Vb, dyt, dyb))) {
    # Calculate alfa for a given system composition
    alfa <- Vt * dyt / (Vt * dyt + Vb * dyb)
  } else{
    AQSys.err("14")
  }
  # the system of equations below was uses the method described by Merchuk
  # in the manuscript Merchuk, J.C., B.A. Andrews, and J.A. Asenjo. (1998),
  # Aqueous two-phase systems for protein separation: Studies on phase inversion.
  # Journal of Chromatography B: Biomedical Sciences and Applications. 711(1-2):
  #p. 285-293. to calculate tieline composition using merchuk's equation.
  # the lines below set the equation's system
  sys <- function(x) {
    F1 <- P1 * exp(P2 * x[2] ^ 0.5 - P3 * x[2] ^ 3) - x[1]
    F2 <- P1 * exp(P2 * x[4] ^ 0.5 - P3 * x[4] ^ 3) - x[3]
    F3 <- (Ym / alfa) - ((1 - alfa) / alfa) * x[3] - x[1]
    F4 <- (Xm / alfa) - ((1 - alfa) / alfa) * x[4] - x[2]
    #
    c(F1 = F1, F2 = F2, F3 = F3, F4 = F4)
  }
  # solve the system of equation for a given set of guess and restricting of positive
  # only results
  (sysres <- multiroot(
    f = sys,
    start = c(1, 0, 0, 1),
    positive = TRUE
  ))
  # Calculate the tieline length and store it in sysres under the TLL alias
  sysres$TLL <- sqrt((sysres$root[1] - sysres$root[3]) ^ 2 + (sysres$root[2] - sysres$root[4])^2)
  # set var name for root results (phase's composition for a given tieline)
  names(sysres$root) <- c("YT", "XT", "YB", "XB")
  # calculate and store tie-line's slope
  sysres$S <- (sysres$root["YT"] - sysres$root["YB"]) / (sysres$root["XT"] - sysres$root["XB"])
  # removing Slope's header to make easier its retrieve
  names(sysres$S) <- NULL
  # return all calculated parameters
  sysres
}
