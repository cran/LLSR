#' @rdname AQSys.CritPoint
#' @name AQSys.CritPoint
#' @title ATPS Critical Point Calculation
#' @description Merchuk et al. described a very straightforward method to calculate the critical composition of a given Binodal
#' curve and its Tielines.
#' @details Using the binodal data, tieline's Slopes (S), the composition of bottom-rich component in the bottom phase (XB)
#' and an equation which stablish a relatioship between them, this function returns the critical composition of the binodal
#' under study. When used within a iterative function, mrchk.tielines() can be used to obtain TLL and S and therefore
#' calculate the critical composition.
#' @method AQSys CritPoint
#' @export AQSys.CritPoint
#' @export
#' @param dataSET - Binodal Experimental data that will be used in the nonlinear fit. [type:data.frame]
#' @param tldata - A data.frame with two columns containing a set of Tieline's Slopes (S)
#' and its bottom-rich component composition in the bottom phase (XB). [type:data.frame]
#' @param method - Binodal Experimental data that will be used in the nonlinear fit. [type:string]
#' "algebraic" - Uses the critical point own definition to set up constraints and solve a system of equations.
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
