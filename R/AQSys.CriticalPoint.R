#' @rdname AQSys.CritPoint
#' @name AQSys.CritPoint
#' @title ATPS Critical Point Calculation
#' @description This function implements methods available in current literature to calculate an ATPS critical point based on its experimental data.
#' @details The Critical Point is one in which both the composition and volume of the phases become equal, and the tie-line length (TLL) tends to zero. 
#' Thus, the methods here implemented the methods decribed by KAUL, A (2000) to calculated a theoretical critical point. 
#' @method AQSys CritPoint
#' @export AQSys.CritPoint
#' @export
#' @param dataSET - Binodal Experimental data that will be used in the nonlinear fit. [type:data.frame]
#' @param tldata - A data.frame with two columns containing a set of Tieline's Slopes (S)
#' and its bottom-rich component composition in the bottom phase (XB). [type:data.frame]
#' @param method - Binodal Experimental data that will be used in the nonlinear fit. [type:string]
#' "algebraic" - Uses the critical point own definition to set up constraints and solve a system of equations. Still in development.
#' "numerical" - A number of tie-lines are calculated successively until TLL is close to zero and concentration of components are numerically equal. A constant slope is assumed.
#' "polynomial" - Calculate the intercept point between the chosen mathematical description and a third order polynomial fitting the tie-lines mid-points
#' @param modelName - Mathematical descriptor that will be used for non-linear fitting. 
#' Use AQSysList() to list the available equations. [type:string]
#' @param slope The method assumes all tielines for a given ATPS are parallel, thus only one slope is required. [type:double]
#' @param NP Number of points used to build the fitted curve. Default is 100. [type:integer]
#' @param xmax Maximum value for the Horizontal axis' value (bottom-rich component). [type:double]
#' @param xlbl Plot's Horizontal axis label. [type:string]
#' @param ylbl Plot's Vertical axis label. [type:string]
#' @param Order Defines how the data is organized in the Worksheet. 
#' Use "xy" whether the first column corresponds to the lower phase fraction and "yx" whether the opposite. [type:string]
#' @param ext - False: Return only XC and YC. True: return an extended output result, 
#' including phase diagram plot and an data.frame including the calculated data. [type:boolean]
#' @param ... Additional optional arguments. None are used at present.
#' @return (XC,YC) - The function returns Tieline's Critical Point Composition
#' @examples
#' \dontrun{
#' AQSys.CritPoint(dataSET, tldata)
#' }
#' @references 
#' KAUL, A. The Phase Diagram. In: HATTI-KAUL, R. (Ed.). Aqueous Two-Phase Systems: Methods and Protocols: Humana Press, v.11, 2000. cap. 2, p.11-21.  (Methods in Biotechnology). ISBN 978-0-89603-541-6.
#' (\href{https://doi.org/10.1385/1-59259-028-4:11}{ACS Publications})
AQSys.CritPoint <- function(dataSET,
                       tldata,
                       method,
                       modelName = "merchuk",
                       slope = NULL,
                       NP = 100,
                       xmax = 30,
                       xlbl = "",
                       ylbl = "",
                       Order = "xy",
                       ext = FALSE,
                       ...) {
  switch(
    method,
    "polynomial" = {
      OUTPUT <- crit_point_poly(
        dataSET = dataSET,
        tldt = tldata,
        modelName = modelName,
        xmax = xmax,
        xlbl = xlbl,
        ylbl = ylbl,
        Order = tolower(Order),
        ext = ext
      )
    },
    "numerical" = {
      OUTPUT <-
        crit_point_seq(
          dataSET = dataSET,
          modelName = modelName,
          NP = NP,
          slope = slope,
          xmax = xmax,
          xlbl = xlbl,
          ylbl = ylbl,
          Order = tolower(Order),
          ext = ext
        )
    },
    "algebraic" = {
      OUTPUT <- crit_point_eqsys(
        dataSET = dataSET,
        tldt = tldata,
        modelName = modelName,
        xmax = xmax,
        xlbl = xlbl,
        ylbl = ylbl,
        Order = tolower(Order),
        ext = ext
      )
    },
    AQSys.err("0")
  )
  cat("\n")
  return(suppressWarnings(OUTPUT))
}
