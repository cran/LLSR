options(digits = 14)
#' @rdname LLSRxy
#' @title Dataset aggregation Function
#' @description This function returns a data.frame from given two unidimensionals vectors, XC and YX.
#' @details The function returns a data.frame after merging two unidimensional vectors.
#' @export LLSRxy
#' @param FirstCol - Component X's concentration in the TOP Y-rich phase.
#' @param SecondCol - Component Y's concentration in the BOTTOM X-rich phase.
#' @param ColDis Defines how the data is organized in the Worksheet. Use "xy" whether the first column corresponds to the lower phase fraction and "yx" whether the opposite.
#' @return XYdt - Standard bidimensional data.frame used in most of functions available in this package.
#' @examples
#' #
#' Xdt<-peg4kslt[,1]
#' #
#' Ydt<-peg4kslt[,1]
#' #
#' LLSRxy(Xdt,Ydt)
LLSRxy <- function(FirstCol, SecondCol, ColDis='xy') {
  # convert and name variables accordingly into vectors
  if (tolower(ColDis) == "xy") {
    xc <- as.vector(as.numeric(sub(",", ".", FirstCol, fixed = TRUE)))
    yc <- as.vector(as.numeric(sub(",", ".", SecondCol, fixed = TRUE)))
  } else{
    xc <- as.vector(as.numeric(sub(",", ".", SecondCol, fixed = TRUE)))
    yc <- as.vector(as.numeric(sub(",", ".", FirstCol, fixed = TRUE)))
  }
  # and combine them into a dataframe
  XYdt <- data.frame(XC = xc, YC = yc)
  # Remove NA's
  XYdt <- XYdt[complete.cases(XYdt), ]
  #
  if ((nrow(XYdt) > 0) && (max(XYdt)<=1)){
    XYdt <- XYdt * 100
  }
  #return data silently - should it be Visible or hidden?
  invisible(XYdt)
}



  
