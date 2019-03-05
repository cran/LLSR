#' @rdname AQSys.LevArmRule
#' @title Lever-Arm Rule - tie-line's Composition Calculation
#' @description Merchuk et al. described a very straightforward method to calculate the concentration of each component in the
#' tieline giving only its global composition and phase's properties (such as volume and density). Here this method is implemented and generalized for multiple mathematical descriptors.
#' @details Using any implemented binodal data mathematical descriptor, the global composition of a chosen tieline and its phases properties.
#' @method AQSys LevArmRule
#' @export AQSys.LevArmRule
#' @export
#' @param dataSET - Binodal Experimental data that will be used in the nonlinear fit
#' @param modelName - Character String specifying the nonlinear empirical equation to fit data.
#' The default method uses Merchuk's equation. Other mathematical descriptors can be listed using AQSysList().
#' @param Xm - Component X's concentration in the tieline's global composition.
#' @param Ym - Component Y's concentration in the tieline's global composition.
#' @param Vt - Tieline's TOP phase volume.
#' @param Vb - Tieline's BOTTOM phase volume.
#' @param dyt - Tieline's TOP phase density
#' @param dyb - Tieline's BOTTOM phase density
#' @param WT - ATPS upper phase weight
#' @param WB - ATPS bottom phase weight
#' @param byW - Use weight (TRUE) or volume and density (FALSE) during lever arm rule calculation.
#' @param Order Defines how the data is organized in the Worksheet. Use "xy" whether the first column corresponds to the lower phase fraction and "yx" whether the opposite.
# ' @param maxiter	- A positive integer specifying the maximum number of iterations allowed.
#' @param ... Additional optional arguments. None are used at present.
#' @return The function returns the Critical Point (X,Y), Tieline Length (TLL), Tieline's Equivolume point (xVRe2o,yVRe2o),
#' and Tieline's Slope.
#' @examples
#' \dontrun{
#' AQSys.LevArmRule(dataSET, Xm, Ym, Vt, Vb, dyt, dyb, WT, WB, byW = FALSE)
#' }
#' @references 
#' MERCHUK, J. C.; ANDREWS, B. A.; ASENJO, J. A. Aqueous two-phase systems for protein separation: Studies on phase inversion. Journal of Chromatography B: Biomedical Sciences and Applications, v. 711, n. 1-2, p. 285-293,  1998. ISSN 0378-4347.
#' (\href{https://www.doi.org/10.1016/s0378-4347(97)00594-x}{ScienceDIrect})
#' 
AQSys.LevArmRule <-
  function(dataSET,
           modelName = "merchuk", 
           Xm,
           Ym,
           Vt = NULL,
           Vb = NULL,
           dyt = NULL,
           dyb = NULL,
           WT = NULL,
           WB = NULL,
           byW = TRUE,
           Order = 'xy',
           ...) {
  #
  dataSET <- toNumeric(dataSET, Order)
  # Fit dataSET data to the chosen equation and subsequently store it in Smmry
  if (modelName %in% names(AQSysList(TRUE))) {
    ans <- do.call(modelName, list(dataSET))
  } else{
    AQSys.err("0")
  }
  Smmry <- summary(ans)
  #
  # extract regression parameters from Smmry
  PARNumber <- AQSysList(TRUE)[[modelName]]
  PARs <- Smmry$coefficients
  #
  if (byW & !is.null(c(WT, WB))) {
    # Calculate alfa for a given system composition
    alfa <- WT / (WT + WB)
  } else{
    AQSys.err("13")
  }
  #
  if (!byW & !is.null(c(Vt, Vb, dyt, dyb))) {
    # Calculate alfa for a given system composition
    alfa <- Vt * dyt / (Vt * dyt + Vb * dyb)
  } else if (!byW){
    AQSys.err("14")
  }
  #
  BnFn <- mathDescPair(modelName)
  FnTerms <- paste(rep("seq(2, 4, 2)", (PARNumber - 1)), collapse = ", ")
  FnChar <- 'sprintf(gsub("\\\\$2", "%d", sprintf(gsub("\\\\$1", "%d", BnFn), seq(1, 3, 2))), $terms)'
  FNs <- eval(parse(text = gsub('\\$terms', FnTerms, FnChar)))
  #
  sys <- function(x) {
    F1 <- eval(parse(text = FNs[1]))
    F2 <- eval(parse(text = FNs[2]))
    F3 <- (Ym / alfa) - ((1 - alfa) / alfa) * x[3] - x[1]
    F4 <- (Xm / alfa) - ((1 - alfa) / alfa) * x[4] - x[2]
    #
    c(F1 = F1, F2 = F2, F3 = F3, F4 = F4)
  }
  #
  sysres <- multiroot(
    f = sys,
    start = c(10, 10, 10, 10),
    positive = TRUE
  )
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
