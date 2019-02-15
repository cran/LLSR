#' @import rootSolve minpack.lm
#' @rdname AQSysOthmer
#' @title Othmer's Equation - Tieline's correlation
#' @description Othmer's equation to correlate tieline's data applying the lever's rule.
#' @param dataSET - Tieline Experimental data that will be used in the nonlinear fit
#' @param Order Defines how the data is organized in the Worksheet. Use "xy" whether the first column corresponds to the lower phase fraction and "yx" whether the opposite.
# ' @param maxiter	- A positive integer specifying the maximum number of iterations allowed.
#' @param ... Additional optional arguments. None are used at present.
#' @export AQSysOthmer
#' @return Parameters A, B and Statistical data
#' @references 
#' OTHMER, D.; TOBIAS, P. Liquid-Liquid Extraction Data - The Line Correlation. Industrial & Engineering Chemistry, v. 34, n. 6, p. 693-696, 1942/06/01 1942. ISSN 0019-7866. 
#' (\href{https://pubs.acs.org/doi/abs/10.1021/ie50390a600}{ACS Publications})
#' @examples
#' # dataSET is a data.frame which contains series of Tieline's mass fraction
#' # (upper-rich component, bottom-rich component and water)
#' # Each column in the data.frame represents a series of one component mass fraction
#' # For example, an empty data.frame for four tielines can be obtaining using:
#' dataSET<-matrix(NA,nrow=4,ncol=6)
#' # Variables order must follows the sequence presented below:
#' # "mfXt","mfYt","mfXb","mfYb","mfWt","mfWb"
#' # In which: mf stands for mass fraction; X and Y for the component
#' # rich in bottom and upper phase, respectively; t or b for top and
#' # bottom phases.
#' # Then you just need to load the data.frame in the function:
#' \dontrun{
#' AQSysOthmer(dataSET)
#'}
AQSysOthmer <- function(dataSET, Order, ...) {
  # store tieline data into a dataframe variable. It might be a better approach check if
  # user stored it in a dataframe and if not trigger an error.
  dataSET <- as.data.frame(dataSET)
  # tieline data is a set of mass fractions of all systems components obtained
  # experimentally for the system's upper and bottom phase.
  # the line bellow set the dataset header
  names(dataSET) <- c("mfXt", "mfYt", "mfXb", "mfYb")
  # Check what range (0-1 or 0-100) is used and standardize
  if ((nrow(dataSET) > 0) && (max(dataSET[1, ]) <= 1)) {
    dataSET <- dataSET * 100
  }
  # the system below will calculate n and K for a given set of tielines
  suppressWarnings(
    FFn <- nlsLM(
      log((100 - mfXb) / mfXb) ~ A + B * log((100 - mfYt) / mfYt),
      start = list(A = 1, B = 1),
      data = dataSET,
      control = nls.lm.control(maxiter = 25)
    )
  )
  # return all calculated parameters
  FFn
}
#' @name AQSysBancroft
#' @title Bancroft's Potential Equation - tie-line's correlation
#' @description Bancroft's equation to correlate tie-line's data.
#' @references 
#' TUBIO, G.  et al. Liquid-liquid equilibrium of the Ucon 50-HB5100/sodium citrate aqueous two-phase systems. Separation and Purification Technology, v. 65, n. 1, p. 3-8,  2009. ISSN 1383-5866. 
#' (\href{https://www.sciencedirect.com/science/article/pii/S1383586608000361}{ScienceDirect})
#' @export AQSysBancroft
#' @param dataSET - Tieline Experimental data that will be used in the nonlinear fit
#' @param Order Defines how the data is organized in the Worksheet. Use "xy" whether the first column corresponds to the lower phase fraction and "yx" whether the opposite.
# ' @param maxiter	- A positive integer specifying the maximum number of iterations allowed.
#' @param ... Additional optional arguments. None are used at present.
#' @return Parameters k1, r and Statistical data
#' @examples
#' # dataSET is a data.frame which contains series of Tieline's mass fraction
#' # (upper-rich component, bottom-rich component and water)
#' # Each column in the data.frame represents a series of one component mass fraction
#' # For example, an empty data.frame for four tielines can be obtaining using:
#' dataSET <- matrix(NA,nrow=4,ncol=6)
#' # Variables order must follows the sequence presented below:
#' # "mfXt","mfYt","mfXb","mfYb","mfWt","mfWb"
#' # In which: mf stands for mass fraction; X and Y for the component
#' # rich in bottom and upper phase, respectively; t or b for top and
#' # bottom phases and W for water.
#' # Then you just need to load the data.frame in the function:
#' \dontrun{
#' AQSysBancroft(dataSET, Order="xy")
#'}
AQSysBancroft <- function(dataSET, Order, ...) {
  # store tie-line data into a dataframe variable. It might be a better approach check if
  # user stored it in a dataframe and if not trigger an error.
  # dataSET <- as.data.frame(dataSET)
  dataSET <- toNumeric(dataSET, Order)
  # tie-line data is a set of mass fractions of all systems components obtained
  # experimentally for the system's upper and bottom phase.
  # the line bellow set the dataset header
  names(dataSET) <- c("mfXt","mfYt","mfXb","mfYb")
  # Check what range (0-1 or 0-100) is used and standardize
  if ((nrow(dataSET) > 0) && (max(dataSET[1, ]) <= 1)) {
    dataSET <- dataSET * 100
  }
  # the system below will calculate r and K1 for a given set of tielines
  suppressWarnings(
    FFn <- nlsLM(
      log(((100 - mfXb - mfYb) / mfXb)) ~ log(k1 * (((100 - mfXt - mfYt) / mfYt
      ) ^ r)),
      start = list(r = 1, k1 = 1),
      data = dataSET,
      control = nls.lm.control(maxiter = 25)
    )
  )
  # return all calculated parameters
  FFn
}
